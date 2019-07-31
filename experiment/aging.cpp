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
#include "common/system.hpp"
#include "common/timer.hpp"
#include "graph/edge.hpp"
#include "graph/edge_stream.hpp"
#include "library/interface.hpp"
#include "configuration.hpp"

using namespace common;
using namespace std;

namespace experiment {

/*****************************************************************************
 *                                                                           *
 * Debug                                                                     *
 *                                                                           *
 *****************************************************************************/
#define DEBUG
#define COUT_DEBUG_FORCE(msg) { LOG("[Aging::" << __FUNCTION__ << "] [" << concurrency::get_thread_id() << "] " << msg); }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif

Aging::Aging(std::shared_ptr<library::UpdateInterface> interface, std::shared_ptr<graph::WeightedEdgeStream> stream, uint64_t num_operations, int64_t num_threads) :
        m_interface(interface), m_stream(stream),
        m_vertices_final(move(* (stream->vertex_table().get()) ) ),
        m_num_operations_total(num_operations), m_num_threads(num_threads), m_num_edges(stream->num_edges()),
        m_max_vertex_id(std::max<uint64_t>(/* avoid overflows */ stream->max_vertex_id(), stream->max_vertex_id() * /* noise */ 2)) {
    LOG("Aging experiment initialised. Number of threads: " << m_num_threads << ", number of operations: " << m_num_operations_total);
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

void Aging::all_workers_execute(const vector<unique_ptr<AgingThread>>& workers, AgingOperation operation){
    vector<future<void>> sync;
    for(auto& w: workers){ sync.push_back( w->execute(operation) ); }
    for(auto& s: sync) { s.get(); }
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
        int64_t part_length = (m_max_vertex_id +1) / num_partitions;
        int64_t odd_partitions = (m_max_vertex_id +1) % num_partitions;
        int64_t start = 0;

        for(int64_t part_id = 0; part_id < num_partitions; part_id++){
            int64_t length = part_length + (part_id < odd_partitions);
            partitions[part_id % m_num_threads].emplace_back(start, length);
            start += length; // next partition
        }
    }

    // start the threads
    LOG("[Aging] Initialising " << m_num_threads << " threads ...");
    vector<thread> threads; threads.reserve(m_num_threads);
    vector<unique_ptr<AgingThread>> workers; workers.reserve(m_num_threads);
    for(int i = 0; i < m_num_threads; i++){
        workers.push_back(make_unique<AgingThread>( this, partitions[i], /* worker_id */ i +1 ));
        threads.emplace_back(&AgingThread::main_thread, workers[i].get());
    }

    // init all threads
    all_workers_execute(workers, AgingOperation::START);
    all_workers_execute(workers, AgingOperation::COMPUTE_FINAL_EDGES);
    m_stream.reset(); // we don't need anymore the list of edges

    // execute the experiment
    LOG("[Aging] Experiment started!");
    Timer timer; timer.start();
    all_workers_execute(workers, AgingOperation::EXECUTE_EXPERIMENT);
    timer.stop();
    LOG("[Aging] Experiment done!");

    // stop all threads
    LOG("[Aging] Waiting for all worker threads to stop...");
    all_workers_execute(workers, AgingOperation::STOP);
    for(auto& t : threads) t.join();
    LOG("[Aging] Worker threads terminated");

    // remove all vertices that are not present in the final graph
    LOG("[Aging] Removing the vertices in surplus, that should not appear in the final graph ...");
    auto lst_vertices = m_vertices_present.lock_table();
    for(auto vertex : lst_vertices){
        if(!m_vertices_final.contains(vertex.first)){ // FIXME to check
            m_interface->delete_vertex(vertex.first);
        }
    }
    lst_vertices.unlock();
    LOG("[Aging] Surplus vertices removed");

    interface->on_thread_destroy(0); // main

    return timer.duration<chrono::microseconds>();
}


/*****************************************************************************
 *                                                                           *
 *  Aging thread                                                             *
 *                                                                           *
 *****************************************************************************/
Aging::AgingThread::AgingThread(Aging* instance, const std::vector<AgingPartition>& partitions, int worker_id) : m_instance(instance), m_interface(instance->m_interface.get()), m_worker_id(worker_id), m_partitions(partitions) {
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

    return future;
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

    assert(m_current_operation == AgingOperation::STOP && "Invalid state, it should still be in the loop");
    m_interface->on_thread_destroy(m_worker_id);
    m_callback.set_value_at_thread_exit();
}

void Aging::AgingThread::main_experiment(){
    // constants
    const uint64_t max_number_edges = static_cast<uint64_t>(m_instance->m_expansion_factor * m_instance->m_num_edges);
    const int64_t num_total_ops = static_cast<int64_t>(m_instance->m_num_operations_total);
//    double prob_insert_final = static_cast<double>(m_edges.size()) / (static_cast<double>(num_total_ops)); // probability of inserting an edge
//    LOG("prob_insert_final: " << prob_insert_final);

    // internal generator
    mt19937_64 random { random_device{}() };
    uniform_real_distribution<double> uniform{ 0., 1. };
    uniform_int_distribution<uint64_t> genrndsrc {0, m_num_src_vertices_in_partitions -1 };
    uniform_int_distribution<uint64_t> genrnddst {0, m_instance->m_max_vertex_id };
    uniform_real_distribution<double> genrndweight{0, configuration().max_weight() };

    int64_t num_ops_done = 0;
    while( (num_ops_done = m_instance->m_num_operations_performed.fetch_add(m_instance->m_granularity) ) < num_total_ops ){

        // shall we perform a burst of insertions or deletions ?
        if(m_edges2remove.empty() /* There are no edges to remove */ ||
                m_interface->num_edges() < max_number_edges /* the size of the current graph is no more than (exp_factor)x of the final graph */){

            // insert `m_granularity' edges then
            for(int64_t i = 0, end = m_instance->m_granularity; i < end; i++){
                double prob_insert_final = static_cast<double>(missing_edges_final()) / (num_total_ops - num_ops_done);
                if ( uniform(random) < prob_insert_final){ // insert from the final graph
                    assert(m_final_edges_current_position < m_edges.size());

                    auto edge = m_edges[m_final_edges_current_position];
                    m_final_edges_current_position++;

                    // if this edge has been previously inserted remove it
                    auto raw_edge = edge.edge();
                    COUT_DEBUG("[Prob: " << prob_insert_final << "] ADD_EDGE FINAL: " << edge);
                    auto res = m_edges_already_inserted.insert_or_assign(raw_edge, /* final */ true);
                    if(! res.second ){ /* the edge was already present */
                        m_interface->delete_edge(raw_edge);
                    }
                    insert_edge(edge);

                } else {
                    // insert a random edge (noise)
                    uint64_t src_id = src_rel2abs(genrndsrc(random));
                    uint64_t dst_id = genrnddst(random);
                    double weight = genrndweight(random);
                    graph::WeightedEdge edge { src_id, dst_id, weight };
                    auto raw_edge = edge.edge();
                    auto res = m_edges_already_inserted.find(raw_edge);
                    bool do_insert = true;
                    bool is_already_registered = false;
                    if(res == m_edges_already_inserted.end()){ // this edge is not already present
                        m_edges_already_inserted[raw_edge] = false;
                    } else if (!res->second){
                        // already present, but it's not final, overwrite its value
                        is_already_registered = true;
                        m_interface->delete_edge(raw_edge);
                    } else {
                        // already present and final, that is it belongs to the final graph
                        do_insert = false;
                    }

                    if(do_insert){
//                        COUT_DEBUG("[Prob: " << prob_insert_final << "] ADD_EDGE TEMP: " << edge);
                        insert_edge(edge);
                        if(!is_already_registered) m_edges2remove.append(raw_edge); /* otherwise already present in the list of edges to remove */
                    } /* else {
                        global_operation_count--;
                    }*/
                }
            }
        } else {
            // perform a burst of deletions
            for(uint64_t i = 0, end = std::min<uint64_t>(m_instance->m_granularity, m_edges2remove.size()); i < end && !m_edges2remove.empty(); i++){
                remove_temporary_edge();
            }
        } // end if (burst of insertions or deletions)

    } // end while (operation count)

    // insert the missing edges from the final graph
    for( ; m_final_edges_current_position < m_edges.size(); m_final_edges_current_position++){
        auto edge = m_edges[m_final_edges_current_position];

        // if this edge has been previously inserted remove it
        auto raw_edge = edge.edge();
        auto res = m_edges_already_inserted.insert_or_assign(raw_edge, /* final */ true);

        COUT_DEBUG("ADD_EDGE FINAL [2]: " << edge);
        if(! res.second ){ /* the edge was already present */
            m_interface->delete_edge(raw_edge);
        }

        insert_edge(edge);
    }

    // remove all edges that do not belong to the final graph
    while(!m_edges2remove.empty()){
        remove_temporary_edge();
    }
}

void Aging::AgingThread::insert_edge(graph::WeightedEdge edge){
    m_instance->insert_edge(edge);
}

void Aging::AgingThread::remove_temporary_edge(){
    assert(!m_edges2remove.empty());

    bool removed = false;
    while(!m_edges2remove.empty() && !removed){
        auto edge = m_edges2remove[0]; m_edges2remove.pop();
        assert(m_edges_already_inserted.find(edge) != m_edges_already_inserted.end() && "This should be in the list of the edges inserted");
        auto it = m_edges_already_inserted.find(edge);
        if(it->second == false /* this is not an edge of the final graph */ ){
//            COUT_DEBUG("DELETE_EDGE TEMP: " << edge);
            m_interface->delete_edge(edge);
            m_edges_already_inserted.erase(it);
            removed = true;
        }
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

bool Aging::AgingThread::vertex_belongs(uint64_t vertex_id) const {
    for(auto& p: m_partitions){
        uint64_t start = p.m_start;
        uint64_t end = start + p.m_length;
        if(vertex_id >= start && vertex_id < end){
            return true;
        }
    }

    return false;
}


} // namespace experiment
