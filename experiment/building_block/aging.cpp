/**
 * Copyright (C) 2019 Dean De Leo, email: dleo[at]cwi.nl
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "aging.hpp"

#include <cassert>
#include <random>
#include <thread>
#include <vector>

#include "common/optimisation.hpp"
#include "common/timer.hpp"
#include "graph/edge.hpp"
#include "graph/edge_stream.hpp"
#include "library/interface.hpp"
#include "configuration.hpp"

using namespace common;
using namespace std;

namespace experiment {

Aging::Aging(std::shared_ptr<library::UpdateInterface> interface) : Aging(interface, make_shared<graph::WeightedEdgeStream>()){ }
Aging::Aging(std::shared_ptr<library::UpdateInterface> interface, std::shared_ptr<graph::WeightedEdgeStream> stream) : Aging(interface, stream, configuration().num_threads(THREADS_WRITE)) { }
Aging::Aging(std::shared_ptr<library::UpdateInterface> interface, std::shared_ptr<graph::WeightedEdgeStream> stream, uint64_t num_operations, int64_t num_threads) : m_interface(interface), m_stream(stream),
        m_num_operations_total(num_operations),
        m_num_threads(num_threads),
        m_vertices_final(move(* (stream->vertex_table().get()) ) ),
        m_vertex_id(std::max<uint64_t>(/* avoid overflows */ stream->max_vertex_id(), stream->max_vertex_id() * /* noise */ 2)) {
    LOG("Aging experiment initialised!");
}

void Aging::set_expansion_factor(double factor){
    if(factor < 1) INVALID_ARGUMENT("the expansion factor must be >= 1, instead the value given is: " << factor);
    LOG("[Aging] Expansion factor set to: " << factor);
    m_expansion_factor = factor;
}

void Aging::set_operation_granularity(uint64_t granularity){
    if(granularity < 1) INVALID_ARGUMENT("the granularity given must be > 0: " << granularity);
    LOG("[Aging] Granularity set to: " << granularity);
    m_granularity = granularity;
}

void Aging::insert_edge(graph::WeightedEdge edge){
    auto interface = m_interface.get();
    if(m_vertices_present.insert(edge.source(), true)) interface->add_vertex(edge.source());
    if(m_vertices_present.insert(edge.destination(), true)) interface->add_vertex(edge.destination());

    // the function returns true if the edge has been inserted. Repeat the loop if it cannot insert the edge as one of
    // the vertices is still being inserted by another thread
    while( ! interface->add_edge(edge) ) { /* nop */ }
}

// run the experiment
std::chrono::microseconds Aging::execute(){
    auto interface = m_interface.get();
    interface->on_main_init(m_num_threads + /* main */ 1);
    interface->on_thread_init(0);

    // compute the set of partitions for each worker
    std::vector<AgingPartition> partitions[m_num_threads];
    {
        int64_t num_partitions = m_num_threads * 4;
        int64_t part_length = (m_vertex_id +1) / num_partitions;
        int64_t odd_partitions = (m_vertex_id +1) % num_partitions;
        int64_t start = 0;

        for(int64_t part_id = 0; part_id < num_partitions; part_id++){
            int64_t length = part_length + (part_id < odd_partitions);
            partitions[part_id % m_num_threads].emplace_back(start, length);
            start += length; // next partition
        }
    }

    // start the threads
    LOG("Initialising " << m_num_threads << " threads ...");
    vector<thread> threads; threads.reserve(m_num_threads);
    vector<AgingThread> workers; workers.reserve(m_num_threads);
    for(int64_t i = 0; i < m_num_threads; i++){
        workers.emplace_back(this, partitions[i], /* worker_id */ i +1);
        threads.emplace_back(&AgingThread::main_thread, workers[i]);
    }

    // init all threads
    vector<future<void>> sync;
    for(auto& w: workers){ sync.push_back( w.execute(AgingOperation::COMPUTE_FINAL_EDGES) ); }
    for(auto& s: sync){ s.get(); }
    sync.clear();
    m_stream.reset(); // we don't need anymore the list of edges

    // execute the experiment
    Timer timer; timer.start();
    for(auto& w: workers){ sync.push_back( w.execute(AgingOperation::EXECUTE_EXPERIMENT) ); }
    LOG("Experiment started!");
    for(auto& s: sync){ s.get(); }
    timer.stop();
    LOG("Experiment done!");
    sync.clear();

    // stop all threads
    LOG("Waiting for all worker threads to stop...");
    for(auto& w: workers){ sync.push_back( w.execute(AgingOperation::STOP) ); }
    for(auto& s: sync){ s.get(); }
    sync.clear();
    for(auto& t : threads) t.join();
    LOG("Worker threads terminated");

    // remove all vertices that are not present in the final graph
    LOG("Removing the vertices in surplus, that should not appear in the final graph ...");
    auto lst_vertices = m_vertices_present.lock_table();
    for(auto vertex : lst_vertices){
        if(m_vertices_final.find(vertex.first) == end(m_vertices_final)){
            m_interface->delete_vertex(vertex.first);
        }
    }
    lst_vertices.unlock();
    LOG("Surplus vertices removed");

    interface->on_thread_destroy(0); // main

    return timer.duration<chrono::microseconds>();
}


/*****************************************************************************
 *                                                                           *
 *  Aging thread                                                             *
 *                                                                           *
 *****************************************************************************/
Aging::AgingThread::AgingThread(Aging* instance, std::vector<AgingPartition> partitions, int worker_id) : m_instance(instance), m_interface(instance->m_interface.get()), m_worker_id(worker_id), m_partitions(partitions) {
    // compute the number of vertices we are responsible to handle
    for(auto& p : m_partitions) { m_num_src_vertices_in_partitions += p.m_length; }
}

Aging::AgingThread::~AgingThread(){

}

std::future<void> Aging::AgingThread::execute(AgingOperation operation){
    auto current_op = m_current_operation;
    if(current_op != AgingOperation::NONE) ERROR("Invalid state: " << (int64_t) operation << ". The worker id " << m_worker_id << " is already busy performing another operation");

    m_callback = promise<void>{ };
    auto future = m_callback.get_future();

    compiler_barrier();
    m_current_operation = operation;
    compiler_barrier();

    return std::move(future);
}

void Aging::AgingThread::main_thread(){
    while(true){
        while(m_current_operation == AgingOperation::NONE) { } // active wait
        if(m_current_operation == AgingOperation::STOP) break; // exit from the loop

        switch(m_current_operation){
        case AgingOperation::START:
            m_interface->on_thread_init(m_worker_id);
            break;
        case AgingOperation::COMPUTE_FINAL_EDGES: {
            m_edges.clear();
            m_final_edges_current_position = 0;
            graph::WeightedEdgeStream* stream = m_instance->m_stream.get();
            for(uint64_t i = 0, end = stream->num_edges(); i < end; i++){
                auto edge = stream->get(i);
                if(vertex_belongs(edge.source())) m_edges.push_back(edge);
            }
        } break;
        case AgingOperation::EXECUTE_EXPERIMENT: {
            main_experiment();
        } break;
        default:
            assert(false && "Invalid operation");
        }

        m_current_operation = AgingOperation::NONE; // move on
        compiler_barrier();
        m_callback.set_value();
    }

    assert(m_current_operation = AgingOperation::STOP && "Invalid state, it should still be in the loop");
    m_interface->on_thread_destroy(m_worker_id);
    m_callback.set_value_at_thread_exit();
}

void Aging::AgingThread::main_experiment(){
    // internal generator
    mt19937_64 random( random_device() );
    uniform_real_distribution<double> uniform{ 0., 1. };
    uniform_int_distribution<uint64_t> genrndsrc {0, m_num_src_vertices_in_partitions };
    uniform_int_distribution<uint64_t> genrnddst {0, m_instance->m_stream->max_vertex_id() };
    uniform_int_distribution<int32_t> genrndweight{0, numeric_limits<int32_t>::max() };

    int64_t global_operation_count = 0;
    while( (global_operation_count = m_instance->m_num_operations_performed.fetch_add(m_instance->m_granularity)) < m_instance->m_num_operations_total ){

        // shall we perform a burst of insertions or deletions ?
        if(m_edges2remove.empty() /* There are no edges to remove */ ||
                m_interface->num_edges() < static_cast<uint64_t>(m_expansion_factor * m_stream->num_edges()) /* the size of the current graph is no more than the 50% of the final graph */){

            // insert `m_granularity' edges then
            for(uint64_t i = 0, end = m_instance->m_granularity; i < end; i++){

                // shall we insert an edge from the final graph or a random edge?
                if ( uniform(random) < ( static_cast<double>(missing_edges_final()) / (missing_edges_final() + (m_instance->m_num_operations_total - global_operation_count) ) ) ) {
                    // insert from the final graph

                    assert(m_final_edges_current_position < m_edges.size());


                    auto edge = m_edges[m_final_edges_current_position];
                    m_final_edges_current_position++;

                    // if this edge has been previously inserted remove it
                    auto raw_edge = edge.edge();
                    auto res = m_edges_already_inserted.insert_or_assign(raw_edge, /* final */ true);
                    if(! res.second /* the edge was already present */){
                        m_interface->delete_edge(raw_edge);
                    }
                    insert_edge(edge);

                } else {
                    // insert a random edge (noise)
                    uint64_t src_id = src_rel2abs(genrndsrc(random));
                    uint64_t dst_id = genrnddst(random);
                    uint32_t weight = static_cast<uint32_t>( genrndweight(random) );
                    graph::WeightedEdge edge { src_id, dst_id, weight };
                    auto raw_edge = edge.edge();
                    auto res = m_edges_already_inserted.find(raw_edge);
                    bool do_insert = true;
                    if(res == m_edges_already_inserted.end()){ // this edge is not already present
                        m_edges_already_inserted[raw_edge] = false;
                    } else if (!res->second){
                        // already present, but it's not final, overwrite its value
                        m_interface->delete_edge(raw_edge);
                    } else {
                        // already present and final, that is it belongs to the final graph
                        do_insert = false;
                    }

                    if(do_insert){
                        insert_edge(edge);
                        m_edges2remove.append(raw_edge);
                    } /* else {
                        global_operation_count--;
                    }*/

                }

                global_operation_count++;
            }

        } else {
            // perform a burst of deletions

            for(uint64_t i = 0, end = std::min<uint64_t>(m_instance->m_granularity, m_edges2remove.size()); i < end; i++){
                remove_temporary_edge();
            }

        } // end if (burst of insertions or deletions)

    } // end while (operation count)


    // remove all edges that do not belong to the final graph
    while(!m_edges2remove.empty()){
        remove_temporary_edge();
    }
}

void Aging::AgingThread::insert(graph::WeightedEdge edge){
    m_instance->insert_edge(edge);
}

void Aging::AgingThread::remove_temporary_edge(){
    assert(!m_edges2remove.empty());
    auto edge = m_edges2remove[0]; m_edges2remove.pop();
    remove_temporary_edge(edge);
}

void Aging::AgingThread::remove_temporary_edge(graph::Edge edge){
    assert(m_edges_already_inserted.find(edge) != m_edges_already_inserted.end() && "This should be in the list of the edges inserted");
    auto it = m_edges_already_inserted.find(edge);
    if(it->second == false /* this is not an edge of the final graph */ ){
        m_interface->delete_edge(edge);
        m_edges_already_inserted.erase(it);
    }
}

int64_t Aging::AgingThread::missing_edges_final() const {
    assert(m_final_edges_current_position <= m_edges.size());
    return static_cast<int64_t>(m_edges.size() - m_final_edges_current_position);
}

uint64_t Aging::AgingThread::src_rel2abs(uint64_t relative_vertex_id) const {
    assert(! m_partitions.empty() && "There are no partitions for this worker");
    int64_t count = static_cast<int64_t>(relative_vertex_id);
    int64_t i = 0, sz = m_partitions.size();
    while(i < sz){
        if(count - m_partitions[i].m_length < 0){
            return m_partitions[i].m_start + count;
        } else {
            count -= m_partitions[i].m_length;
            i++;
        }
    }

    assert(false && "Invalid vertex ID");
    ERROR("Invalid relative vertex ID: " << relative_vertex_id);
}


} // namespace experiment
