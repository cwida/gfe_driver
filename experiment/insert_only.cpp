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

#include "insert_only.hpp"

#include <atomic>
#include <chrono>
#include <future>
#include <thread>
#include <vector>

#include "common/database.hpp"
#include "common/system.hpp"
#include "common/timer.hpp"
#include "configuration.hpp"
#include "library/interface.hpp"
#include "third-party/libcuckoo/cuckoohash_map.hh"

using namespace common;
using namespace std;

namespace experiment {

InsertOnly::InsertOnly(std::shared_ptr<library::UpdateInterface> interface, std::shared_ptr<graph::WeightedEdgeStream> graph, int64_t num_threads, bool measure_latency) :
        m_interface(interface), m_stream(graph), m_num_threads(num_threads), m_measure_latency(measure_latency) {
    if(m_num_threads == 0) ERROR("Invalid number of threads: " << m_num_threads);
}

bool InsertOnly::is_static_scheduler() const {
    return m_schedule_chunks == 0;
}

void InsertOnly::set_round_robin_scheduler(uint64_t granularity) {
    if(granularity == 0) INVALID_ARGUMENT("The granularity of a chunk size cannot be zero");
    m_schedule_chunks = granularity;
}

void InsertOnly::set_static_scheduler(){
    m_schedule_chunks = 0;
}

void InsertOnly::set_batch_size(uint64_t batch_size){
    if(m_measure_latency) ERROR("Batch size not supported when measuring the latency");
    m_batch_size = batch_size;
}

// Execute an update at the time
static void run_one_by_one(library::UpdateInterface* interface, graph::WeightedEdgeStream* graph, cuckoohash_map<uint64_t, bool>& vertices, uint64_t start, uint64_t end){
    for(uint64_t pos = start; pos < end; pos++){
        auto edge = graph->get(pos);

        if(vertices.insert(edge.source(), true)) interface->add_vertex(edge.source());
        if(vertices.insert(edge.destination(), true)) interface->add_vertex(edge.destination());

        // the function returns true if the edge has been inserted. Repeat the loop if it cannot insert the edge as one of
        // the vertices is still being inserted by another thread
        while( ! interface->add_edge(edge) ) { /* nop */ } ;
    }
}

// Execute an update at the time, measure the latency of each insertion
static void run_one_by_one(library::UpdateInterface* interface, graph::WeightedEdgeStream* graph, cuckoohash_map<uint64_t, bool>& vertices, uint64_t start, uint64_t end, uint64_t* __restrict latencies){
    for(uint64_t pos = start; pos < end; pos++){
        auto edge = graph->get(pos);

        if(vertices.insert(edge.source(), true)) interface->add_vertex(edge.source());
        if(vertices.insert(edge.destination(), true)) interface->add_vertex(edge.destination());

        // the function returns true if the edge has been inserted. Repeat the loop if it cannot insert the edge as one of
        // the vertices is still being inserted by another thread
        auto t0 = chrono::steady_clock::now();
        while( ! interface->add_edge(edge) ) {
            t0 = chrono::steady_clock::now();
        };


        common::compiler_barrier();
        auto t1 = chrono::steady_clock::now();
        latencies[pos] = chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
    }
}

// Send batches of updates
static void run_batched(library::UpdateInterface* interface, graph::WeightedEdgeStream* graph, cuckoohash_map<uint64_t, bool>& vertices, uint64_t start, uint64_t end, uint64_t batch_size){
    unique_ptr<library::UpdateInterface::SingleUpdate[]> ptr_batch { new library::UpdateInterface::SingleUpdate[batch_size] };
    library::UpdateInterface::SingleUpdate* __restrict batch = ptr_batch.get();
    uint64_t batch_index = 0;

    for(uint64_t pos = start; pos < end; pos++){
        auto edge = graph->get(pos);

        if(vertices.insert(edge.source(), true)) interface->add_vertex(edge.source());
        if(vertices.insert(edge.destination(), true)) interface->add_vertex(edge.destination());

        if(batch_index == batch_size){
            interface->batch(batch, batch_index);
            batch_index = 0; // reset
        }

        library::UpdateInterface::SingleUpdate* __restrict update = batch + batch_index;
        update->m_source = edge.m_source;
        update->m_destination = edge.m_destination;
        assert(edge.m_weight >= 0 && "Weights can only be non-negative");
        update->m_weight = edge.m_weight;
        batch_index++;
    }

    // flush
    interface->batch(batch, batch_index);
}

// this wasn't well thought, I know
static void run_sequential(library::UpdateInterface* interface, graph::WeightedEdgeStream* graph, cuckoohash_map<uint64_t, bool>& vertices, uint64_t start, uint64_t end, uint64_t batch_size, uint64_t* latencies = nullptr){
    if(batch_size <= 0){ // send one update at the time
        if(latencies == nullptr){ // do not measure the latency of each insertion
            run_one_by_one(interface, graph, vertices, start, end);
        } else { // measure the latency of each insertion
            run_one_by_one(interface, graph, vertices, start, end, latencies);
        }
    } else { // send updates in batches
        assert(latencies == nullptr && "Cannot measure the latency of insertions in batch mode");
        run_batched(interface, graph, vertices, start, end, batch_size);
    }
}

void InsertOnly::execute_static(void* cb, uint64_t* latencies){
    auto graph = m_stream.get();
    auto& vertices = *( reinterpret_cast<cuckoohash_map<uint64_t, bool>*>(cb) );


    int64_t edges_per_threads = graph->num_edges() / m_num_threads;
    int64_t odd_threads = graph->num_edges() % m_num_threads;
    int64_t start = 0;

    vector<thread> threads;

    for(int64_t i = 0; i < m_num_threads; i++){
        int64_t length = edges_per_threads + (i < odd_threads);

        threads.emplace_back([this, graph, &vertices, latencies](int64_t start, int64_t length, int thread_id){
            concurrency::set_thread_name("Worker #" + to_string(thread_id));

            auto interface = m_interface.get();
            interface->on_thread_init(thread_id);
            run_sequential(interface, graph, vertices, start, start + length, m_batch_size, latencies /* do not shift */);
            interface->on_thread_destroy(thread_id);

        }, start, length, static_cast<int>(i));

        start += length;
    }

    // wait for all threads to complete
    for(auto& t : threads) t.join();
}

void InsertOnly::execute_round_robin(void* cb, uint64_t* latencies){
    auto& vertices = *( reinterpret_cast<cuckoohash_map<uint64_t, bool>*>(cb) );
    ASSERT(m_schedule_chunks > 0 && "Chunk size == 0");

    vector<thread> threads;

    atomic<uint64_t> start_chunk_next = 0;

    for(int64_t i = 0; i < m_num_threads; i++){
        threads.emplace_back([this, &vertices, &start_chunk_next, latencies](int thread_id){
            concurrency::set_thread_name("Worker #" + to_string(thread_id));

            auto interface = m_interface.get();
            auto graph = m_stream.get();
            uint64_t start;
            const uint64_t size = graph->num_edges();

            interface->on_thread_init(thread_id);

            while( (start = start_chunk_next.fetch_add(m_schedule_chunks)) < size ){
                uint64_t end = std::min<uint64_t>(start + m_schedule_chunks, size);
                run_sequential(interface, graph, vertices, start, end, m_batch_size, latencies /* do not shift */);
            }

            interface->on_thread_destroy(thread_id);

        }, static_cast<int>(i));
    }

    // wait for all threads to complete
    for(auto& t : threads) t.join();
}

chrono::microseconds InsertOnly::execute() {
    // check which vertices have been already inserted
    cuckoohash_map<uint64_t, bool> vertices;
    vertices.max_num_worker_threads(thread::hardware_concurrency());
    std::unique_ptr<uint64_t[]> ptr_latencies { m_measure_latency ? new uint64_t[m_stream->num_edges()] : nullptr };
    uint64_t* latencies = ptr_latencies.get();

    m_interface->on_main_init(m_num_threads);

    Timer timer;
    timer.start();

    if(is_static_scheduler()){
        execute_static(&vertices, latencies);
    } else {
        execute_round_robin(&vertices, latencies);
    }
    timer.stop();
    LOG("Insertions performed with " << m_num_threads << " threads in " << timer);
    m_time_insert = timer.microseconds();

    m_interface->on_thread_init(0);
    timer.start();
    m_interface->build();
    timer.stop();
    m_interface->on_thread_destroy(0);
    m_time_build = timer.microseconds();
    if(m_time_build > 0){
        LOG("Build time: " << timer);
    }

    m_interface->on_main_destroy();

    LOG("Edge stream size: " << m_stream->num_edges() << ", num edges stored in the graph: " << m_interface->num_edges() << ", match: " << (m_stream->num_edges() == m_interface->num_edges() ? "yes" : "no"));

    if(m_measure_latency){
        LOG("Computing the measured latencies...");
        m_latencies = details::LatencyStatistics::compute_statistics(latencies, m_stream->num_edges());
        LOG("Measured latencies: " << m_latencies);
    }

    return chrono::microseconds{ m_time_insert + m_time_build };
}

void InsertOnly::save() {
    assert(configuration().db() != nullptr);
    auto db = configuration().db()->add("insert_only");
    db.add("scheduler", is_static_scheduler() ? "static" : "round_robin");
    db.add("insertion_time", m_time_insert); // microseconds
    db.add("build_time", m_time_build); // microseconds
    db.add("num_edges", m_stream->num_edges());

    if(m_measure_latency) m_latencies.save("inserts");
}

} // namespace experiment
