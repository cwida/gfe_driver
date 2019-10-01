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


#include "aging2_worker.hpp"

#include <cassert>
#include <chrono>
#include <iostream>
#include <mutex>
#include <random>
#include <thread>
#include <utility>

#include "common/error.hpp"
#include "common/system.hpp"
#include "common/timer.hpp"
#include "experiment/aging2_experiment.hpp"
#include "graph/edge_stream.hpp"
#include "library/interface.hpp"
#include "aging2_master.hpp"
#include "configuration.hpp"

using namespace common;
using namespace std;

/*****************************************************************************
 *                                                                           *
 * Debug                                                                     *
 *                                                                           *
 *****************************************************************************/
extern mutex _log_mutex [[maybe_unused]];
//#define DEBUG
#define COUT_DEBUG_FORCE(msg) { scoped_lock<mutex> lock(_log_mutex); cout << "[Aging2Worker::" << __FUNCTION__ << "] [" << concurrency::get_thread_id() << ", worker_id: " << m_worker_id << "] " << msg << endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif

/*****************************************************************************
 *                                                                           *
 * Init                                                                      *
 *                                                                           *
 *****************************************************************************/
namespace experiment::details {

Aging2Worker::Aging2Worker(Aging2Master& master, int worker_id) : m_master(master), m_library(m_master.parameters().m_library.get()), m_worker_id(worker_id), m_task(Task::IDLE) {
    assert(m_library != nullptr);

    // compute the range of this partition for this worker
    int64_t length = m_master.parameters().m_stream->num_edges() / m_master.parameters().m_num_threads;
    int64_t modulo = m_master.parameters().m_stream->num_edges() % m_master.parameters().m_num_threads;
    m_partition_start = std::min<int64_t>(worker_id, modulo) * (length +1) +
            std::max<int64_t>(0, static_cast<int64_t>(worker_id) - modulo) * length;
    if(worker_id < modulo) length ++;
    m_partition_end = m_partition_start + length;

    // start the background thread
    start();
}

Aging2Worker::~Aging2Worker(){
    stop();
}

void Aging2Worker::start(){
   assert(m_thread.joinable() == false && "The background thread should not be already running at this point");
   assert(m_task == Task::IDLE);

   m_task = Task::START;
   unique_lock<mutex> lock(m_mutex);
   m_thread = std::thread(&Aging2Worker::main_thread, this);
   m_condvar.wait(lock, [this](){return m_task == Task::IDLE; });

   assert(m_thread.joinable() == true);
}

void Aging2Worker::stop(){
    assert(m_thread.joinable() == true && "Already stopped");

    unique_lock<mutex> lock(m_mutex);
    assert(m_task == Task::IDLE && "Service busy");
    m_task = Task::STOP;
    lock.unlock();
    m_condvar.notify_all();
    m_thread.join();
}

/*****************************************************************************
 *                                                                           *
 * Interface                                                                 *
 *                                                                           *
 *****************************************************************************/

void Aging2Worker::execute_updates(){
    set_task_async(Task::EXECUTE_UPDATES);
}

void Aging2Worker::remove_artificial_vertices(){
    set_task_async(Task::REMOVE_ARTIFICIAL_VERTICES);
}

void Aging2Worker::set_task_async(Task task){
    { // restrict the scope
        scoped_lock<mutex> lock(m_mutex);
        assert(m_task == Task::IDLE && "Service busy");
        if(m_task != Task::IDLE) ERROR("Background thread already busy performing another operation");
        m_task = task;
    }
    m_condvar.notify_all();
}

void Aging2Worker::wait(){
    unique_lock<mutex> lock(m_mutex);
    m_condvar.wait(lock, [this](){ return m_task == Task::IDLE; });
}

/*****************************************************************************
 *                                                                           *
 * Background thread                                                         *
 *                                                                           *
 *****************************************************************************/

void Aging2Worker::main_thread(){
    COUT_DEBUG("Worker started, partition: [" << m_partition_start << ", " << m_partition_end << "), size = " << (m_partition_end - m_partition_start));

    m_library->on_thread_init(m_worker_id);

    bool terminate = false;
    Task task; // current task
    do {
        { // fetch the next task to execute
            unique_lock<mutex> lock(m_mutex);
            assert(m_task != Task::IDLE && "Incorrect state");
            m_task = Task::IDLE;
            m_condvar.notify_all();
            // wait for the next task to execute
            m_condvar.wait(lock, [this](){ return m_task != Task::IDLE; });
            task = m_task;
        }

        switch(task){
        case Task::IDLE:
            assert(0 && "Invalid operation");
            break;
        case Task::START:
            assert(0 && "This operation is reserved only for starting the service, it should not occur anymore at this point");
            break;
        case Task::STOP:
            terminate = true;
            break;
        case Task::EXECUTE_UPDATES:
            main_execute_updates();
            break;
        case Task::REMOVE_ARTIFICIAL_VERTICES:
            main_remove_artificial_vertices();
            break;
        }
    } while (!terminate);

    m_library->on_thread_destroy(m_worker_id);

    // not really necessary, only present for consistency ..
    unique_lock<mutex> lock(m_mutex);
    assert(m_task != Task::IDLE && "Incorrect state");
    m_task = Task::IDLE;
    m_condvar.notify_all();

    // we're done
    COUT_DEBUG("Worker terminated");
}

void Aging2Worker::main_execute_updates(){
    const int64_t num_total_ops = m_master.num_operations_total();
    const uint64_t max_num_edges = m_master.max_num_edges(); // max num edges that the graph can reach in the intermediate graph
    const bool report_progress = m_master.parameters().m_report_progress;
    uniform_int_distribution<uint64_t> genrndsrc {0, m_master.distr_max_value() -1 }; // incl.
    uniform_real_distribution<double> genrndweight{0, m_master.parameters().m_max_weight }; // incl.

    // heuristics to bump up the probability of inserting a final edge due to multiple threads and deletions
    const double prob_bump = 1.0 * m_master.parameters().m_num_threads;

    int64_t num_ops_done = 0;
    int lastset_coeff = 0;

    while( (num_ops_done = m_master.m_num_operations_performed.fetch_add( granularity() )) < num_total_ops ){

        // Report progress
        if(report_progress && static_cast<int>(100.0 * num_ops_done/num_total_ops) > m_master.m_last_progress_reported){
            m_master.m_last_progress_reported = 100.0 * num_ops_done/num_total_ops;
#if defined(DEBUG)
            LOG("[thread: " << common::concurrency::get_thread_id() << ", worker_id: " << m_worker_id << "] "
                    "Progress: " << num_ops_done << "/" << num_total_ops << " (" << 100.0 * num_ops_done/num_total_ops << "%), "
                    "edges final graph: " <<  m_partition_offset << "/" << (m_partition_end - m_partition_start) << " (" << (100.0 * m_partition_offset/ (m_partition_end - m_partition_start)) << " %), "
                    "hash load factor: " << (100.0 * m_master.m_edges_present.load_factor()) << " %"
            );
#else // just report the percentages
            LOG("[thread: " << common::concurrency::get_thread_id() << ", worker_id: " << m_worker_id << "] "
                    "Progress: " << 100.0 * num_ops_done/num_total_ops << "%, "
                    "edges final graph: " <<  (100.0 * m_partition_offset/ (m_partition_end - m_partition_start)) << " %"
            );
#endif
        }

        // shall we perform a burst of insertions or deletions ?
        if(m_edges2remove.empty() /* There are no edges to remove */ ||
                m_library->num_edges() < max_num_edges /* the size of the current graph is no more than (exp_factor)x of the final graph */){

            // insert `m_granularity' edges then
            for(uint64_t i = 0, end = granularity(); i < end; i++){
                double prob_insert_final = prob_bump * static_cast<double>(missing_edges_final()) / (num_total_ops - num_ops_done);
                if ( m_uniform(m_random) < prob_insert_final ){ // insert from the final graph
                    assert(m_partition_start + m_partition_offset < m_partition_end && "No edges to insert from the final graph");

                    auto edge = m_master.parameters().m_stream->get(m_partition_start + m_partition_offset);
                    if(! m_master.is_directed() && edge.source() > edge.destination() ) edge.swap_src_dst();

                    // if this edge has been previously inserted remove it
                    auto raw_edge = edge.edge();
                    if( m_master.m_edges_present.uprase_fn(raw_edge, [this, edge](bool& value){
                        // the edge is already present in the graph
                        assert(value == false && "This must be an artificial edge");
                        graph_remove_edge(edge.edge(), true);
                        graph_insert_edge(edge);
                        value = true; // the edge now belongs to the final graph
                        return false; // false => do not erase, but keep the entry in the hash map
                    },  true ) ){
                        // this edge was not present in the graph
                        graph_insert_edge(edge);
                    }

                    assert(m_master.m_edges_present.contains(raw_edge));
                    assert(m_master.m_edges_present.find(raw_edge) == /* this is the value associated to the key */ true);

                    m_partition_offset++;
                } else {
                    // insert a random edge (noise)
                    bool inserted = false;

                    do {
                        // generate the source vertex
                        uint64_t source_rank = m_master.distr_vertex_offset( genrndsrc(m_random) );
                        uint64_t source_vertex_id = m_master.m_vertices_freq [ source_rank ].m_vertex_id;
                        uint64_t source_Cfreq0 = source_rank > 0 ? m_master.m_vertices_freq [ source_rank -1 ].m_frequency : 0;
                        uint64_t source_Cfreq1 = m_master.m_vertices_freq [ source_rank ].m_frequency;
                        uint64_t source_freq = source_Cfreq1 - source_Cfreq0;

                        // generate the destination vertex
                        uint64_t destination_random_sample = uniform_int_distribution<uint64_t> {0, m_master.distr_max_value() -1 - source_freq }(m_random);
                        if(destination_random_sample >= source_Cfreq0) destination_random_sample += source_freq; // avoid generating an edge src -> src
                        uint64_t destination_rank = m_master.distr_vertex_offset(destination_random_sample);
                        assert(source_rank != destination_rank);
                        uint64_t destination_vertex_id = m_master.m_vertices_freq[ destination_rank ].m_vertex_id;
                        COUT_DEBUG("source: " << source_vertex_id << ", destination: " << destination_vertex_id);
                        assert(source_vertex_id != destination_vertex_id);

                        double weight = genrndweight(m_random);
                        graph::WeightedEdge edge { source_vertex_id, destination_vertex_id, weight };
                        if(edge.source() > edge.destination()) edge.swap_src_dst();
                        auto raw_edge = edge.edge();

                        if(m_master.m_edges_present.insert(raw_edge, /* edge belongs to final graph ? */ false)){
                            graph_insert_edge(edge);
                            m_edges2remove.append(raw_edge); /* otherwise already present in the list of edges to remove */
                            inserted = true;
                        }

                    } while (!inserted);

                }

            }
        } else {
            // perform a burst of deletions
            for(uint64_t i = 0, end = std::min<uint64_t>(granularity(), m_edges2remove.size()); i < end && !m_edges2remove.empty(); i++){
                graph_remove_temporary_edge();
            }
        } // end if (burst of insertions or deletions)

        // report how long it took to perform 1x, 2x, ... updates w.r.t. to the size of the final graph
        int aging_coeff = (num_ops_done + granularity()) / m_master.parameters().m_stream->num_edges();
        if(aging_coeff > lastset_coeff){
            if( m_master.m_last_time_reported.compare_exchange_strong(/* updates lastset_coeff */ lastset_coeff, aging_coeff) ){
                uint64_t duration = chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - m_master.m_time_start ).count();
                m_master.m_reported_times[aging_coeff -1] = duration;
            }
        }

    } // end while (operation count)

    // insert the missing edges from the final graph
    COUT_DEBUG("Processed edges: " << m_partition_offset << " / " << (m_partition_end - m_partition_start) << " (" << (100.0 * m_partition_offset/(m_partition_end - m_partition_start)) << " %)");
    for( ; m_partition_start + m_partition_offset < m_partition_end; m_partition_offset++){
        auto edge = m_master.parameters().m_stream->get(m_partition_start + m_partition_offset);

        // if this edge has been previously inserted remove it
        auto raw_edge = edge.edge();
        if( m_master.m_edges_present.uprase_fn(raw_edge, [this, edge](bool& value){
            // the edge is already present in the graph
            assert(value == false && "This must be an artificial edge");
            graph_remove_edge(edge.edge(), /* force remove */ true);
            graph_insert_edge(edge);
            value = true; // the edge now belongs to the final graph
            return false; // false => do not erase, but keep the entry in the hash map
        },  true ) ){
            // this edge was not present in the graph
            graph_insert_edge(edge);
        }
    }

    // remove all edges that do not belong to the final graph
    while(!m_edges2remove.empty()){
        graph_remove_temporary_edge();
    }
}

void Aging2Worker::main_remove_artificial_vertices(){
    assert(m_master.m_vertices2remove != nullptr && "Array not initialised in the master node");

    uint64_t* __restrict vertices2remove = m_master.m_vertices2remove;
    uint64_t vertices2remove_sz = m_master.m_results.num_artificial_vertices();

    // compute the range of this partition for this worker
    int64_t length = vertices2remove_sz / m_master.parameters().m_num_threads;
    int64_t modulo = vertices2remove_sz % m_master.parameters().m_num_threads;
    int64_t start = std::min<int64_t>(m_worker_id, modulo) * (length +1) +
            std::max<int64_t>(0, static_cast<int64_t>(m_worker_id) - modulo) * length;
    if(m_worker_id < modulo) length ++;

    for(uint64_t i = start, end = start + length; i < end; i++){
        m_library->remove_vertex(vertices2remove[i]);
    }
}

/*****************************************************************************
 *                                                                           *
 * Utility functions                                                         *
 *                                                                           *
 *****************************************************************************/

void Aging2Worker::graph_insert_vertex(uint64_t vertex_id){
    if(m_master.m_vertices_present.insert(vertex_id, true)){
        COUT_DEBUG("insert vertex: " << vertex_id);
        m_library->add_vertex(vertex_id);
    }
}

void Aging2Worker::graph_insert_edge(graph::WeightedEdge edge){
    // be sure that the vertices source & destination are already present
    graph_insert_vertex(edge.source());
    graph_insert_vertex(edge.destination());

    if(!m_master.is_directed() && m_uniform(m_random) < 0.5) edge.swap_src_dst(); // noise
    COUT_DEBUG("edge: " << edge);

    // the function returns true if the edge has been inserted. Repeat the loop if it cannot insert the edge as one of
    // the vertices is still being inserted by another thread
    while ( ! m_library->add_edge(edge) ) { /* nop */ };
}

void Aging2Worker::graph_remove_edge(graph::Edge edge, bool force){
    if(!m_master.is_directed() && m_uniform(m_random) < 0.5) edge.swap_src_dst(); // noise
    COUT_DEBUG("edge: " << edge);
    if(!force){
        m_library->remove_edge(edge);
    } else { // force = true
        while( ! m_library->remove_edge(edge) ) /* nop */ ;
    }
}

void Aging2Worker::graph_remove_temporary_edge(){
    assert(!m_edges2remove.empty());

    bool removed = false;
    while(!m_edges2remove.empty() && !removed){
        auto edge = m_edges2remove[0]; m_edges2remove.pop();

        m_master.m_edges_present.erase_fn(edge, [this, &removed, edge](bool& is_final_graph){
            if(!is_final_graph){
                graph_remove_edge(edge);
                removed = true;
                return true; // true => remove the element from the hash table
            } else {
                return false; // false => the edge belongs to the final graph, do not remove it
            }
        });
    }
}

uint64_t Aging2Worker::granularity() const{
    return m_master.parameters().m_worker_granularity;
}

int64_t Aging2Worker::missing_edges_final() const {
    assert(m_partition_start + m_partition_offset <= m_partition_end);
    return m_partition_end - (m_partition_start + m_partition_offset);
}

} // namespace
