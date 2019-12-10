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
#include "details/build_thread.hpp"
#include "configuration.hpp"
#include "library/interface.hpp"
#include "third-party/libcuckoo/cuckoohash_map.hh"

using namespace common;
using namespace gfe::experiment::details;
using namespace std;

namespace gfe::experiment {

InsertOnly::InsertOnly(std::shared_ptr<gfe::library::UpdateInterface> interface, std::shared_ptr<gfe::graph::WeightedEdgeStream> graph, int64_t num_threads, bool measure_latency) :
        m_interface(interface), m_stream(graph), m_num_threads(num_threads), m_measure_latency(measure_latency) {
    if(m_num_threads == 0) ERROR("Invalid number of threads: " << m_num_threads);
}

void InsertOnly::set_scheduler_granularity(uint64_t granularity) {
    if(granularity == 0) INVALID_ARGUMENT("The granularity of a chunk size cannot be zero");
    m_scheduler_granularity = granularity;
}

void InsertOnly::set_build_frequency(std::chrono::milliseconds millisecs){
    m_build_frequency = millisecs;
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
static void run_one_by_one_with_latency(library::UpdateInterface* interface, graph::WeightedEdgeStream* graph, cuckoohash_map<uint64_t, bool>& vertices, uint64_t start, uint64_t end, uint64_t* __restrict latencies){
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

// this wasn't well thought, I know
static void run_sequential(library::UpdateInterface* interface, graph::WeightedEdgeStream* graph, cuckoohash_map<uint64_t, bool>& vertices, uint64_t start, uint64_t end, uint64_t* latencies = nullptr){
    if(latencies == nullptr){ // do not measure the latency of each insertion
        run_one_by_one(interface, graph, vertices, start, end);
    } else { // measure the latency of each insertion
        run_one_by_one_with_latency(interface, graph, vertices, start, end, latencies);
    }
}

void InsertOnly::execute_round_robin(void* cb, uint64_t* latencies){
    auto& vertices = *( reinterpret_cast<cuckoohash_map<uint64_t, bool>*>(cb) );

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

            while( (start = start_chunk_next.fetch_add(m_scheduler_granularity)) < size ){
                uint64_t end = std::min<uint64_t>(start + m_scheduler_granularity, size);
                run_sequential(interface, graph, vertices, start, end, latencies /* do not shift */);
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

    // re-adjust the scheduler granularity if there are too few insertions to perform
    if(m_stream->num_edges() / m_num_threads < m_scheduler_granularity){
        m_scheduler_granularity = m_stream->num_edges() / m_num_threads;
        if(m_scheduler_granularity == 0) m_scheduler_granularity = 1; // corner case
        LOG("InsertOnly: reset the scheduler granularity to " << m_scheduler_granularity << " edge insertions per thread");
    }

    // Execute the insertions
    m_interface->on_main_init(m_num_threads /* build thread */ +1);
    Timer timer;
    timer.start();
    BuildThread build_service { m_interface , static_cast<int>(m_num_threads), m_build_frequency };
    execute_round_robin(&vertices, latencies);
    build_service.stop();
    timer.stop();
    LOG("Insertions performed with " << m_num_threads << " threads in " << timer);
    m_time_insert = timer.microseconds();
    m_num_build_invocations = build_service.num_invocations();

    // A final invocation of the method #build()
    m_interface->on_thread_init(0);
    timer.start();
    m_interface->build();
    timer.stop();
    m_interface->on_thread_destroy(0);
    m_time_build = timer.microseconds();
    if(m_time_build > 0){
        LOG("Build time: " << timer);
    }
    m_num_build_invocations++;

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
    db.add("scheduler", "round_robin"); // backwards compatibility
    db.add("scheduler_granularity", m_scheduler_granularity); // the number of insertions performed by each thread
    db.add("insertion_time", m_time_insert); // microseconds
    db.add("build_time", m_time_build); // microseconds
    db.add("num_edges", m_stream->num_edges());
    db.add("num_snapshots_created", m_interface->num_levels());
    db.add("num_build_invocations", m_num_build_invocations);
    // missing revision: until 25/Nov/2019
    // version 20191125: build thread, build frequency taken into account, scheduler set to round_robin, removed batch updates
    // version 20191210: difference between num_build_invocations (explicit invocations to #build()) and num_snapshots_created (actual number of deltas created by the impl)
    db.add("revision", "20191210");

    if(m_measure_latency) m_latencies.save("inserts");
}

} // namespace
