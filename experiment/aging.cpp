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
#include <cmath>
#include <condition_variable>
#include <future>
#include <mutex>
#include <random>
#include <thread>
#include <vector>

#include <gmpxx.h> // libgmp

#include "common/database.hpp"
#include "common/optimisation.hpp"
#include "common/system.hpp"
#include "common/timer.hpp"
#include "details/aging_operation.hpp"
#include "details/aging_partition.hpp"
#include "details/aging_thread.hpp"
#include "details/async_batch.hpp"
#include "graph/edge.hpp"
#include "graph/edge_stream.hpp"
#include "graph/vertex_list.hpp"
#include "library/interface.hpp"
#include "configuration.hpp"

using namespace common;
using namespace experiment::details;
using namespace std;

/*****************************************************************************
 *                                                                           *
 * Debug                                                                     *
 *                                                                           *
 *****************************************************************************/
extern mutex _log_mutex [[maybe_unused]];
#define DEBUG
#define COUT_DEBUG_FORCE(msg) { scoped_lock<mutex> lock(_log_mutex); cout << "[Aging::" << __FUNCTION__ << "] [" << concurrency::get_thread_id() << "] " << msg << endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif

/*****************************************************************************
 *                                                                           *
 * Aging                                                                     *
 *                                                                           *
 *****************************************************************************/
namespace experiment {

Aging::Aging(std::shared_ptr<library::UpdateInterface> interface, std::shared_ptr<graph::WeightedEdgeStream> stream, double mult_num_operations, int64_t num_threads) :
        Aging(interface,stream, stream->num_edges() * mult_num_operations, num_threads, interface->is_directed(), stream->max_weight()) { }

Aging::Aging(std::shared_ptr<library::UpdateInterface> interface, std::shared_ptr<graph::WeightedEdgeStream> stream, uint64_t num_operations, int64_t num_threads, bool graph_is_directed, double max_weight) :
        m_interface(interface), m_stream(stream),
        m_vertices_final(move(* (stream->vertex_table().get()) ) ),
        m_num_operations_total(num_operations), m_num_threads(num_threads), m_num_edges(stream->num_edges()), m_num_vertices(m_vertices_final.size()),
        m_max_vertex_id_final(stream->max_vertex_id()),
        m_is_directed(graph_is_directed), m_max_weight(max_weight){
    LOG("Aging experiment initialised. Number of threads: " << m_num_threads << ", number of operations: " << m_num_operations_total);
    m_max_vertex_id_artificial = m_num_vertices * m_ef_vertices;
#if defined(DEBUG)
    set_report_progress(true);
#endif

    // 1024 is a hack to avoid issues with small graphs
    m_reported_times = new uint64_t[ std::max<uint64_t>( ::ceil (static_cast<double>(num_operations) / m_num_edges ) + 1, 1024 ) ]();
}

Aging::~Aging(){
    delete[] m_reported_times; m_reported_times = nullptr;
}

void Aging::set_expansion_factor_edges(double factor){
    if(factor < 1) INVALID_ARGUMENT("the expansion factor must be >= 1, instead the value given is: " << factor);
    LOG("[Aging] Expansion factor for the edges set to: " << factor);
    m_ef_edges = factor;
}

void Aging::set_expansion_factor_vertices(double factor){
    if(factor < 1) INVALID_ARGUMENT("the expansion factor must be >= 1, instead the value given is: " << factor);
    LOG("[Aging] Expansion factor for the vertices set to: " << factor);
    m_ef_vertices = factor;
    m_max_vertex_id_artificial = m_ef_vertices * m_num_vertices;
}

void Aging::set_operation_granularity(uint64_t granularity){
    if(granularity < 1) INVALID_ARGUMENT("the granularity given must be > 0: " << granularity);
    LOG("[Aging] Granularity set to: " << granularity);
    m_granularity = granularity;
}

void Aging::set_batch_size(uint64_t size){
    m_batch_size = size;
}

void Aging::set_report_progress(bool value){
    m_report_progress = value;
}

void Aging::compute_partitions_add_missing(vector<vector<AgingPartition>>* out_partitions) const {
    assert(out_partitions != nullptr);
    if(m_max_vertex_id_final <= m_max_vertex_id_artificial) return; // nop

    int64_t num_partitions = m_num_threads * 8;
    int64_t start = m_max_vertex_id_artificial +1;
    int64_t part_length = (m_max_vertex_id_final - m_max_vertex_id_artificial) / num_partitions;
    int64_t odd_partitions = (m_max_vertex_id_final - m_max_vertex_id_artificial) % num_partitions;
    if(part_length == 0)
        num_partitions = odd_partitions;
    for(int64_t part_id = 0; part_id < num_partitions; part_id++){
        int64_t length = part_length + (part_id < odd_partitions);
        COUT_DEBUG("partition[" << part_id << "] missing for interval: [" << start << ", " << start + length << ")");
        (*out_partitions)[part_id % m_num_threads].emplace_back(start, length);
        start += length; // next partition
    }
}

vector<vector<AgingPartition>> Aging::compute_partitions_directed() const{
    vector<vector<AgingPartition>> partitions ( m_num_threads );
    int64_t num_partitions = m_num_threads * 8;
    int64_t start = 0;

    // partitions for the artificial vertices
    int64_t part_length = (m_max_vertex_id_artificial +1) / num_partitions;
    int64_t odd_partitions = (m_max_vertex_id_artificial +1) % num_partitions;
    for(int64_t part_id = 0; part_id < num_partitions; part_id++){
        int64_t length = part_length + (part_id < odd_partitions);
        partitions[part_id % m_num_threads].emplace_back(start, length);
        start += length; // next partition
    }

    compute_partitions_add_missing(&partitions);
    return partitions;
}

vector<vector<AgingPartition>> Aging::compute_partitions_undirected() const{
    vector<vector<AgingPartition>> partitions ( m_num_threads );
    int64_t num_partitions = m_num_threads * 8;
    mpz_class max_num_edges = mpz_class{m_max_vertex_id_artificial} * mpz_class{m_max_vertex_id_artificial} / 2; // (n-1)*(n-2) /2;
    mpz_class e = max_num_edges / num_partitions; // edges per partition

    LOG("Computing the list of partitions...");
    COUT_DEBUG("num_partitions: " << num_partitions << ", max_vertex_id: " << m_max_vertex_id_artificial << ", max_num_edges: " << max_num_edges << ", edges_per_partition: " << e);

    uint64_t vertex_from = 0;
    int64_t part_id = 0;
    mpz_t tmp_x; mpz_init(tmp_x); // temporary
    constexpr size_t buffer_sz = 1024; char buffer[buffer_sz];

    while(part_id < num_partitions -1 && vertex_from < (m_max_vertex_id_artificial -1)){
        // how many edges can we create from vertex_from ?
        mpz_class S = mpz_class{m_max_vertex_id_artificial +1 -vertex_from} * mpz_class{m_max_vertex_id_artificial -vertex_from} /2;

        // solve the inequation S - [(x^2 +x) /2] < e
        // that is x^2 +x +2(e -s) > 0 => [ -1 + sqrt( 1 - 8(e-s) ) ] / 2
//        uint64_t x = ceil( (-1.0 + sqrt(1ll - 8ll*(e-S)) ) / 2.0 );
        mpz_class sqrt_arg = ( S - e ) * 8 +1;
        mpz_sqrt(tmp_x, sqrt_arg.get_mpz_t());
        mpz_sub_ui(tmp_x, tmp_x, 1u);
        mpz_cdiv_q_ui(tmp_x, tmp_x, 2u); // compute the quotient & round up (the `c' in cdiv)

        // convert the result back to uint64_t
        assert(mpz_sizeinbase(tmp_x, 10) <= buffer_sz && "Conversion overflow");
        mpz_get_str(buffer, 10, tmp_x);
        uint64_t x = strtoull(buffer, nullptr, 10);

        assert(x <= m_max_vertex_id_artificial);
        uint64_t vertex_upto = m_max_vertex_id_artificial - x;
        if(vertex_from == vertex_upto) vertex_upto++; // corner case
        else if(vertex_upto >= m_max_vertex_id_artificial) vertex_upto = (m_max_vertex_id_artificial -1); // corner case

        COUT_DEBUG("partition[" << part_id << "] interval: [" << vertex_from << ", " << vertex_upto << ")");
        partitions[part_id % m_num_threads].emplace_back(vertex_from, vertex_upto - vertex_from);

        vertex_from = vertex_upto; // next iteration
        part_id++;
    }
    mpz_clear(tmp_x);

    COUT_DEBUG("partition[" << part_id << "] interval: [" << vertex_from << ", " << m_max_vertex_id_artificial << "]" );
    partitions[part_id % m_num_threads].emplace_back(vertex_from, m_max_vertex_id_artificial - vertex_from +1);

    compute_partitions_add_missing(&partitions);

    return partitions;
}


bool Aging::is_directed() const {
    return m_is_directed;
}

void Aging::insert_vertex(uint64_t vertex_id){
    if(m_vertices_present.insert(vertex_id, true)){
//        COUT_DEBUG("insert vertex: " << vertex_id);
        m_interface->add_vertex(vertex_id);
    }
}

void Aging::all_workers_execute(const vector<unique_ptr<AgingThread>>& workers, AgingOperation operation){
    vector<future<void>> sync;
    for(auto& w: workers){ sync.push_back( w->execute(operation) ); }
    for(auto& s: sync) { s.get(); }
}

void Aging::log_num_vtx_edges(){
    scoped_lock<mutex> lock(_log_mutex);
    cout << "[Aging] Number of stored vertices: " << m_num_vertices_final_graph << " [match: ";
    if(m_num_vertices == m_num_vertices_final_graph){ cout << "yes"; } else {
        cout << "no, expected " << m_num_vertices;
    }
    cout << "], number of stored edges: " << m_num_edges_final_graph << " [match: ";
    if(m_num_edges == m_num_edges_final_graph){ cout << "yes"; } else {
        cout << "no, expected " << m_num_edges;
    }
    cout << "]" << endl;
}

// run the experiment
std::chrono::microseconds Aging::execute(){
    auto interface = m_interface.get();
    int num_workers =  m_num_threads * 2 + /* main */ 1; // * 2 due to batch async
    interface->on_main_init(num_workers);
    interface->on_thread_init(0);

    // compute the set of partitions for each worker
    vector<vector<AgingPartition>> partitions = is_directed() ? compute_partitions_directed() : compute_partitions_undirected();
    assert(partitions.size() == m_num_threads && "Expected one partition set per thread");

    // start the threads
    LOG("[Aging] Initialising " << m_num_threads << " threads ...");
    vector<thread> threads; threads.reserve(m_num_threads);
    vector<unique_ptr<AgingThread>> workers; workers.reserve(m_num_threads);
    for(int i = 0; i < m_num_threads; i++){
        int worker_id = i * 2 +1; // worker_id is for the AgingThread, and worker_id +1 is for the AsyncBatch
        workers.push_back(make_unique<AgingThread>( this, partitions[i], worker_id ));
        threads.emplace_back(&AgingThread::main_thread, workers[i].get());
    }

    // init all threads
    all_workers_execute(workers, AgingOperation::START);
    all_workers_execute(workers, AgingOperation::COMPUTE_FINAL_EDGES);
    m_stream.reset(); // we don't need anymore the list of edges

    // execute the experiment
    LOG("[Aging] Experiment started!");
    m_last_progress_reported = 0;
    m_last_time_reported = 0; m_time_start = chrono::steady_clock::now();
    Timer timer; timer.start();
    all_workers_execute(workers, AgingOperation::EXECUTE_EXPERIMENT);
    timer.stop();
    LOG("[Aging] Experiment done!");
    LOG("[Aging] Updates performed with " << m_num_threads << " threads in " << timer);
    m_completion_time = timer.microseconds();
    auto completion_time_microsecs = timer.duration<chrono::microseconds>();

    {
        LOG("[Aging] Removing the vertices in excess, that should not appear in the final graph ...");

        timer.start();

        all_workers_execute(workers, AgingOperation::INTERNAL_CLEANUP); // free up some memory
        // remove all vertices that are not present in the final graph
        auto lst_vertices = m_vertices_present.lock_table();
        m_vertices2remove = new uint64_t[lst_vertices.size()];
        uint64_t vertices2remove_sz = 0;
        for(auto vertex : lst_vertices){
            if(!m_vertices_final.contains(vertex.first)){
                m_vertices2remove[vertices2remove_sz] = vertex.first;
                vertices2remove_sz++;
    //            m_interface->remove_vertex(vertex.first);
            }
        }
        lst_vertices.unlock();
        uint64_t items_per_part = vertices2remove_sz / workers.size();
        uint64_t odd_items = vertices2remove_sz % workers.size();
        uint64_t start = 0;
        for(uint64_t i = 0; i < workers.size(); i++){
            uint64_t length = items_per_part + (i < odd_items);
            workers[i]->set_partition_vertices_to_remove(start, length);
            start += length; // next iteration
        }
        all_workers_execute(workers, AgingOperation::REMOVE_VERTICES);
        delete[] m_vertices2remove; m_vertices2remove = nullptr;
        m_num_artificial_vertices = vertices2remove_sz;
        LOG("[Aging] Number of extra vertices: " << m_num_artificial_vertices << ", "
                "expansion factor: " << static_cast<double>(m_num_artificial_vertices + m_num_vertices) / m_num_vertices);
        timer.stop();
        LOG("[Aging] Artificial vertices removed in " << timer);

        m_num_vertices_final_graph = interface->num_vertices();
        m_num_edges_final_graph = interface->num_edges();
        log_num_vtx_edges();
    }


    // stop all threads
    LOG("[Aging] Waiting for all worker threads to stop...");
    all_workers_execute(workers, AgingOperation::STOP);
    for(auto& t : threads) t.join();
    LOG("[Aging] Worker threads terminated");


    interface->on_thread_destroy(0); // main

    return completion_time_microsecs;
}


void Aging::save() {
    assert(configuration().db() != nullptr);
    auto db = configuration().db()->add("aging");
    db.add("granularity", m_granularity);
    db.add("num_threads", m_num_threads);
    db.add("num_updates_requested", m_num_operations_total);
    db.add("num_updates_executed", m_num_operations_performed);
    db.add("num_artificial_vertices", m_num_artificial_vertices);
    db.add("num_vertices_load", m_num_vertices);
    db.add("num_vertices_final", m_num_vertices_final_graph);
    db.add("num_edges_load", m_num_edges);
    db.add("num_edges_final", m_num_edges_final_graph);
    db.add("completion_time", m_completion_time); // microseconds

    LOG("m_last_time_reported: " << m_last_time_reported);
    for(int i = 0, sz = m_last_time_reported; i < sz; i++){
        if(m_reported_times[i] == 0) continue; // missing??
        auto db = configuration().db()->add("aging_intermediate_throughput");
        db.add("aging_coeff", (int64_t) i +1); // 1, 2, 3...
        db.add("completion_time", m_reported_times[i]); // microseconds
    }
}


} // namespace experiment
