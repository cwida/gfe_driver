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

using namespace common;
using namespace gfe::experiment::details;
using namespace std;

namespace gfe::experiment {

InsertOnly::InsertOnly(std::shared_ptr<gfe::library::UpdateInterface> interface, std::shared_ptr<gfe::graph::WeightedEdgeStream> graph, int64_t num_threads) :
        m_interface(interface), m_stream(graph), m_num_threads(num_threads) {
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
static void run_sequential(library::UpdateInterface* interface, graph::WeightedEdgeStream* graph, uint64_t start, uint64_t end){
    for(uint64_t pos = start; pos < end; pos++){
        auto edge = graph->get(pos);
        [[maybe_unused]] bool result = interface->add_edge_v2(edge);
        assert(result == true && "Edge not inserted");
    }
}

void InsertOnly::execute_round_robin(){
    vector<thread> threads;

    atomic<uint64_t> start_chunk_next = 0;

    for(int64_t i = 0; i < m_num_threads; i++){
        threads.emplace_back([this, &start_chunk_next](int thread_id){
            concurrency::set_thread_name("Worker #" + to_string(thread_id));

            auto interface = m_interface.get();
            auto graph = m_stream.get();
            uint64_t start;
            const uint64_t size = graph->num_edges();

            interface->on_thread_init(thread_id);

            while( (start = start_chunk_next.fetch_add(m_scheduler_granularity)) < size ){
                uint64_t end = std::min<uint64_t>(start + m_scheduler_granularity, size);
                run_sequential(interface, graph, start, end);
            }

            interface->on_thread_destroy(thread_id);

        }, static_cast<int>(i));
    }

    // wait for all threads to complete
    for(auto& t : threads) t.join();
}

chrono::microseconds InsertOnly::execute() {
    // re-adjust the scheduler granularity if there are too few insertions to perform
    if(m_stream->num_edges() / m_num_threads < m_scheduler_granularity){
        m_scheduler_granularity = m_stream->num_edges() / m_num_threads;
        if(m_scheduler_granularity == 0) m_scheduler_granularity = 1; // corner case
        LOG("InsertOnly: reset the scheduler granularity to " << m_scheduler_granularity << " edge insertions per thread");
    }

    // Execute the insertions
    m_interface->on_main_init(m_num_threads /* build thread */ +1);
    m_interface->updates_start();
    Timer timer;
    timer.start();
    BuildThread build_service { m_interface , static_cast<int>(m_num_threads), m_build_frequency };
    execute_round_robin();
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

    m_time_build = timer.microseconds();
    if(m_time_build > 0){
        LOG("Build time: " << timer);
    }
    m_num_build_invocations++;

    m_interface->updates_stop();

    LOG("Edge stream size: " << m_stream->num_edges() << ", num edges stored in the graph: " << m_interface->num_edges() << ", match: " << (m_stream->num_edges() == m_interface->num_edges() ? "yes" : "no"));

    m_interface->on_thread_destroy(0);
    m_interface->on_main_destroy();

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
    // version 20200625: rely on #add_edge_v2 to implicitly create the vertices. This should alleviate the footprint of the driver for non scalable implementations
    db.add("revision", "20200625");
}

} // namespace
