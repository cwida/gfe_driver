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
namespace gfe::experiment::details {

Aging2Worker::Aging2Worker(Aging2Master& master, int worker_id) : m_master(master), m_library(m_master.parameters().m_library.get()), m_worker_id(worker_id), m_task{ TaskOp::IDLE, nullptr, 0 } {
    assert(m_library != nullptr);

    // start the background thread
    start();
}

Aging2Worker::~Aging2Worker(){
    stop();

    while(!m_updates.empty()){
        delete m_updates[0]; m_updates.pop();
    }
}

void Aging2Worker::start(){
   assert(m_thread.joinable() == false && "The background thread should not be already running at this point");
   assert(m_task.m_type == TaskOp::IDLE);

   m_task.m_type = TaskOp::START;
   unique_lock<mutex> lock(m_mutex);
   m_thread = std::thread(&Aging2Worker::main_thread, this);
   m_condvar.wait(lock, [this](){return m_task.m_type == TaskOp::IDLE; });

   assert(m_thread.joinable() == true);
}

void Aging2Worker::stop(){
    assert(m_thread.joinable() == true && "Already stopped");

    unique_lock<mutex> lock(m_mutex);
    assert(m_task.m_type == TaskOp::IDLE && "Service busy");
    m_task.m_type = TaskOp::STOP;
    lock.unlock();
    m_condvar.notify_all();
    m_thread.join();
}

/*****************************************************************************
 *                                                                           *
 * Interface                                                                 *
 *                                                                           *
 *****************************************************************************/

void Aging2Worker::load_edges(uint64_t* edges, uint64_t num_edges){
    set_task_async(TaskOp::LOAD_EDGES, edges, num_edges);
}

void Aging2Worker::execute_updates(){
    set_task_async(TaskOp::EXECUTE_UPDATES);
}

void Aging2Worker::remove_vertices(uint64_t* vertices, uint64_t num_vertices){
    set_task_async(TaskOp::REMOVE_VERTICES, vertices, num_vertices);
}

void Aging2Worker::set_latencies(uint64_t* array_insertions, uint64_t* array_deletions){
    set_task_async(TaskOp::SET_ARRAY_LATENCIES, array_insertions, /* hack */ reinterpret_cast<uint64_t>(array_deletions));
}

void Aging2Worker::set_task_async(TaskOp type, uint64_t* payload, uint64_t payload_sz){
    { // restrict the scope
        scoped_lock<mutex> lock(m_mutex);
        assert(m_task.m_type == TaskOp::IDLE && "Service busy");
        if(m_task.m_type != TaskOp::IDLE) ERROR("Background thread already busy performing another operation");
        m_task = Task { type, payload, payload_sz };
    }
    m_condvar.notify_all();
}

void Aging2Worker::wait(){
    unique_lock<mutex> lock(m_mutex);
    m_condvar.wait(lock, [this](){ return m_task.m_type == TaskOp::IDLE; });
}

bool Aging2Worker::wait(const chrono::time_point<chrono::steady_clock>& tp){
    auto now = chrono::steady_clock::now();

    unique_lock<mutex> lock(m_mutex);
    while(m_task.m_type != TaskOp::IDLE && /* timeout ? */ !((tp <= now) || (tp - now) < 1ms) ) {
        m_condvar.wait_until(lock, tp, [this](){ return m_task.m_type == TaskOp::IDLE; });
        now = chrono::steady_clock::now(); // next iteration
    };

    return m_task.m_type == TaskOp::IDLE;
}

/*****************************************************************************
 *                                                                           *
 * Background thread                                                         *
 *                                                                           *
 *****************************************************************************/

void Aging2Worker::main_thread(){
    COUT_DEBUG("Worker started");
    concurrency::set_thread_name("Worker #" + to_string(m_worker_id));

    m_library->on_thread_init(m_worker_id);

    bool terminate = false;
    Task task; // current task
    do {
        { // fetch the next task to execute
            unique_lock<mutex> lock(m_mutex);
            assert(m_task.m_type != TaskOp::IDLE && "Incorrect state");
            m_task = Task{ TaskOp::IDLE, nullptr, 0 };
            m_condvar.notify_all();
            // wait for the next task to execute
            m_condvar.wait(lock, [this](){ return m_task.m_type != TaskOp::IDLE; });
            task = m_task;
        }

        switch(task.m_type){
        case TaskOp::IDLE:
            assert(0 && "Invalid operation");
            break;
        case TaskOp::START:
            assert(0 && "This operation is reserved only for starting the service, it should not occur anymore at this point");
            break;
        case TaskOp::STOP:
            terminate = true;
            break;
        case TaskOp::SET_ARRAY_LATENCIES:
            m_latency_insertions = task.m_payload;
            m_latency_deletions = reinterpret_cast<uint64_t*>(task.m_payload_sz); // hack
            break;
        case TaskOp::LOAD_EDGES:
            main_load_edges(task.m_payload, task.m_payload_sz);
            break;
        case TaskOp::EXECUTE_UPDATES:
            main_execute_updates();
            break;
        case TaskOp::REMOVE_VERTICES:
            main_remove_vertices(task.m_payload, task.m_payload_sz);
            break;
        }
    } while (!terminate);

    m_library->on_thread_destroy(m_worker_id);

    // not really necessary, only present for consistency ..
    unique_lock<mutex> lock(m_mutex);
    assert(m_task.m_type != TaskOp::IDLE && "Incorrect state");
    m_task = Task{ TaskOp::IDLE, nullptr, 0 };
    m_condvar.notify_all();

    // we're done
    COUT_DEBUG("Worker terminated");
}

void Aging2Worker::main_execute_updates(){
    // compute the amount of space used by the vectors in m_updates
    for(uint64_t i = 0; i < m_updates.size(); i++){ m_updates_mem_usage += m_updates[i]->size() * sizeof(gfe::graph::WeightedEdge); } // either size() or capacity()
    COUT_DEBUG("Initial memory footprint: " << m_updates_mem_usage << " bytes");

    const int64_t num_total_ops = m_master.num_operations_total();
    const bool report_progress = m_master.parameters().m_report_progress;
    const bool release_memory = m_master.parameters().m_release_driver_memory;
    // reports_per_ops only affects how often a report is saved in the db, not the report to the stdout
    const double reports_per_ops = m_master.parameters().m_num_reports_per_operations;
    int lastset_coeff = 0;

    for(uint64_t i = 0, end = m_updates.size(); i < end; i++){
        // if we're release the driver's memory, always fetch the first. Otherwise follow the index.
        vector<graph::WeightedEdge>* operations = m_updates[release_memory ? 0 : i];

        uint64_t num_loops = (operations->size() / granularity()) + (operations->size() % granularity() != 0);
        uint64_t start = 0;
        for(uint64_t j = 0; j < num_loops; j++){
            uint64_t end = std::min( start + granularity(), operations->size() );

            // execute a chunk of updates
            graph_execute_batch_updates(operations->data() + start, end - start);

            uint64_t num_ops_done = m_master.m_num_operations_performed.fetch_add( end - start );

            // report progress
            if(report_progress && static_cast<int>(100.0 * num_ops_done/num_total_ops) > m_master.m_last_progress_reported){
                m_master.m_last_progress_reported = 100.0 * num_ops_done/num_total_ops;
                if(!m_master.m_stop_experiment){
                    LOG("[thread: " << ::common::concurrency::get_thread_id() << ", worker_id: " << m_worker_id << "] Progress: " << static_cast<int>(100.0 * num_ops_done/num_total_ops) << "%");
                }
            }

            // report how long it took to perform 1x, 2x, ... updates w.r.t. to the size of the final graph
            int aging_coeff = (static_cast<double>(num_ops_done) / m_master.num_edges_final_graph()) * reports_per_ops;
            if(aging_coeff > lastset_coeff){
                if( m_master.m_last_time_reported.compare_exchange_strong(/* updates lastset_coeff */ lastset_coeff, aging_coeff) ){
                    uint64_t duration = chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - m_master.m_time_start ).count();
                    m_master.m_reported_times[aging_coeff -1] = duration;
                }
            }

            // next iteration
            start = end;
        }

        if(release_memory){
            COUT_DEBUG("Releasing a buffer of cardinality " << operations->size() << ", " << m_updates.size() -1 << " buffers left");
            m_updates_mem_usage -= m_updates[0]->capacity() * sizeof(gfe::graph::WeightedEdge); // update the memory footprint of this worker
            COUT_DEBUG("Memory footprint: " << m_updates_mem_usage << " bytes");
            delete m_updates[0];
            m_updates.pop();
        }
    }
}


void Aging2Worker::main_load_edges(uint64_t* edges, uint64_t num_edges){
    if(m_updates.empty()){ m_updates.append(new vector<graph::WeightedEdge>()); }
    vector<graph::WeightedEdge>* last = m_updates[m_updates.size() -1];

    constexpr uint64_t last_max_sz = (1ull << 22); // 4M
    const uint64_t modulo = m_master.parameters().m_num_threads;

    uint64_t* __restrict sources = edges;
    uint64_t* __restrict destinations = sources + num_edges;
    double* __restrict weights = reinterpret_cast<double*>(destinations + num_edges);

    uniform_real_distribution<double> rndweight{0, m_master.parameters().m_max_weight}; // in [0, max_weight)

    for(uint64_t i = 0; i < num_edges; i++){
        if(static_cast<int>((sources[i] + destinations[i]) % modulo) == m_worker_id){
            if(last->size() > last_max_sz){
                last = new vector<graph::WeightedEdge>();
                m_updates.append(last);
            }

            // counters
            double weight = weights[i];
            if(weight >= 0){
                m_num_edge_insertions++;
            } else {
                m_num_edge_deletions++;
            }

            // generate a random weight
            if(weight == 0.0){
                weight = rndweight(m_random); // in [0, max_weight)
                if(weight == 0.0) weight = m_master.parameters().m_max_weight; // in (0, max_weight]
            }


            last->emplace_back(sources[i], destinations[i], weight);
        }
    }
}

void Aging2Worker::main_remove_vertices(uint64_t* vertices, uint64_t num_vertices){
    auto pardegree = m_master.parameters().m_num_threads;

    for(uint64_t i = 0; i < num_vertices; i++){
        if(static_cast<int>(vertices[i] % pardegree) == m_worker_id){
            COUT_DEBUG("Remove vertex: " << vertices[i]);
            m_library->remove_vertex(vertices[i]);
        }
    }
}

/*****************************************************************************
 *                                                                           *
 * Utility functions                                                         *
 *                                                                           *
 *****************************************************************************/

void Aging2Worker::graph_execute_batch_updates(graph::WeightedEdge* __restrict updates, uint64_t num_updates){
    if(m_latency_insertions == nullptr){
        assert(m_master.parameters().m_measure_latency == false);
        assert(m_latency_deletions == nullptr);
        graph_execute_batch_updates0</* measure latency ? */ false>(updates, num_updates);
    } else {
        assert(m_master.parameters().m_measure_latency == true);
        assert(m_latency_deletions != nullptr);

        graph_execute_batch_updates0</* measure latency ? */ true>(updates, num_updates);
    }
}

template<bool with_latency>
void Aging2Worker::graph_execute_batch_updates0(graph::WeightedEdge* __restrict updates, uint64_t num_updates){
    for(uint64_t i = 0; i < num_updates; i++){
        if(m_master.m_stop_experiment) break; // timeout, we're done

        if(updates[i].m_weight >= 0){ // insertion
            graph_insert_edge<with_latency>(updates[i]);
        } else { // deletion
            graph_remove_edge<with_latency>(updates[i].edge());
        }

        m_num_operations++;
    }
}

template<bool with_latency>
void Aging2Worker::graph_insert_edge(graph::WeightedEdge edge){
    if(!m_master.is_directed() && m_uniform(m_random) < 0.5) edge.swap_src_dst(); // noise
    COUT_DEBUG("edge: " << edge);

    if(with_latency == false){
        // the function returns true if the edge has been inserted. Repeat the loop if it cannot insert the edge as one of
        // the vertices is still being inserted by another thread
        while ( ! m_library->add_edge_v2(edge) ) { /* nop */ };

    } else { // measure the latency of the insertion
        chrono::steady_clock::time_point t0, t1;
        do {
            t0 = chrono::steady_clock::now();
        } while ( ! m_library->add_edge_v2(edge) );
        t1 = chrono::steady_clock::now();

        m_latency_insertions[0] = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
        m_latency_insertions++;
    }

}

template<bool with_latency>
void Aging2Worker::graph_remove_edge(graph::Edge edge, bool force){
    if(!m_master.is_directed() && m_uniform(m_random) < 0.5) edge.swap_src_dst(); // noise
    COUT_DEBUG("edge: " << edge);

    if(with_latency == false){

        if(!force){
            m_library->remove_edge(edge);
        } else { // force = true
            while( ! m_library->remove_edge(edge) ) /* nop */ ;
        }

    } else { // measure the latency of the deletion
        chrono::steady_clock::time_point t0, t1;

        t0 = chrono::steady_clock::now();
        if(!force){
            m_library->remove_edge(edge);
        } else { // force = true
            while( ! m_library->remove_edge(edge) ) /* nop */;
        }
        t1 = chrono::steady_clock::now();

        m_latency_deletions[0] = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
        m_latency_deletions++;
    }

}

uint64_t Aging2Worker::granularity() const{
    return m_master.parameters().m_worker_granularity;
}

uint64_t Aging2Worker::num_insertions() const{
    return m_num_edge_insertions;
}

uint64_t Aging2Worker::num_deletions() const {
    return m_num_edge_deletions;
}

uint64_t Aging2Worker::num_operations() const {
    return m_num_operations;
}

uint64_t Aging2Worker::memory_footprint() const {
    return m_updates_mem_usage;
}

} // namespace
