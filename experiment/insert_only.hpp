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

#include "graph/edge.hpp"
#include "graph/edge_stream.hpp"
#include "library/interface.hpp"

namespace gfe::experiment {

/**
 * This experiment simulates only insertions (and no deletions) in a given graph system. It assumes that the system
 * to evaluate has already been created, and it is empty, that is, it should not already contain any vertices or edges.
 * The experiment inserts all vertices and edges, stored in a given edge stream, in parallel, into the system, and
 * measures the overall throughput.
 */
class InsertOnly {
    std::shared_ptr<gfe::library::UpdateInterface> m_interface; // the library where vertices and edges will be inserted
    std::shared_ptr<gfe::graph::WeightedEdgeStream> m_stream; // the graph to insert
    const int64_t m_num_threads; // the number of threads to use
    std::chrono::milliseconds m_build_frequency {0}; // Continuously create a new snapshot each `m_build_frequency' millisecs (0 = feature disabled)
    uint64_t m_scheduler_granularity = 1ull << 20; // if >0, granularity for the scheduler
    uint64_t m_time_insert = 0; // the amount of time to insert all elements in the database, in microseconds
    uint64_t m_time_build = 0; // the amount of time to build the last snapshot/delta/level in the library, in microseconds
    uint64_t m_num_build_invocations = 0; // number of times the method #build() has been invoked

    // Execute the experiment with the round robin scheduler
    void execute_round_robin();

public:
    // Initialise the experiment
    // @param interface the system to evaluate, already instantiated
    // @param stream the list of edges to insert in the system, possibly already permuted
    // @param num_threads the parallelism degree, that is the number of threads to employ to insert concurrently all edges from the stream
    // @param measure_latency if requested, report the median/min/max/percentiles measured latencies in insertions. Measuring the latency
    //        could introduce additional overhead and decrease the measure for the overall throughput
    InsertOnly(std::shared_ptr<gfe::library::UpdateInterface> interface, std::shared_ptr<gfe::graph::WeightedEdgeStream> stream, int64_t num_threads);

    // Advanced and/or internal parameter, set the granularity of chunks sent by the internal scheduler to the worker threads. The granularity
    // here is given by the number of edge insertions in each chunk.
    void set_scheduler_granularity(uint64_t granularity);

    // Set how frequently create a new snapshot/delta in the library (0 = do not create new snapshots)
    void set_build_frequency(std::chrono::milliseconds millisecs);

    // Execute the experiment
    std::chrono::microseconds execute();

    // Store the results into the database
    void save();
};

} // namespace experiment


