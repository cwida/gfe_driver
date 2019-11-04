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

#pragma once

#include <chrono>
#include <cinttypes>
#include <memory>

#include "details/latency.hpp"
#include "graph/edge.hpp"
#include "graph/edge_stream.hpp"
#include "library/interface.hpp"

namespace experiment {

class InsertOnly {
    std::shared_ptr<library::UpdateInterface> m_interface; // the library where vertices and edges will be inserted
    std::shared_ptr<graph::WeightedEdgeStream> m_stream; // the graph to insert
    const int64_t m_num_threads; // the number of threads to use
    const bool m_measure_latency; // whether to measure and report the latency of each insertion
    uint64_t m_schedule_chunks = 0; // if >0, schedule the edges to insert in round robin fashion, in chunks of the given size
    uint64_t m_time_insert = 0; // the amount of time to insert all elements in the database, in microseconds
    uint64_t m_time_build = 0; // the amount of time to build the snapshot/delta/level in the library, in microseconds
    uint64_t m_batch_size = 0; // send the updates in batches
    details::LatencyStatistics m_latencies; // the latencies measured in the experiment

    // Execute the experiment with the static scheduler
    void execute_static(void* /* opaque */ cb, uint64_t* latencies);

    // Execute the experiment with the round robin scheduler
    void execute_round_robin(void* /* opaque */ cb, uint64_t* latencies);

    // Send the the updates one by one

public:
    InsertOnly(std::shared_ptr<library::UpdateInterface> interface, std::shared_ptr<graph::WeightedEdgeStream> stream, int64_t num_threads, bool measure_latency = false);

    // Check whether to use a round robin or a static scheduler
    bool is_static_scheduler() const;

    // Whether to schedule the edges to the different threads in chunks, with the given granularity of a chunk size
    void set_round_robin_scheduler(uint64_t granularity);

    // Static scheduler, each thread inserts a fixed amount of edges
    void set_static_scheduler();

    // Request to send edge updates in batches of the given size. Vertex insertions will continue to be sent one at the time
    void set_batch_size(uint64_t size);

    // Execute the experiment
    std::chrono::microseconds execute();

    // Store the results into the database
    void save();
};

} // namespace experiment


