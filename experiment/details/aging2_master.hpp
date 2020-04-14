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

#include <memory>
#include <vector>

#include "common/static_index.hpp"
#include "experiment/aging2_result.hpp"
#include "graph/edge.hpp"
#include "third-party/libcuckoo/cuckoohash_map.hh"

// forward declarations
namespace gfe::experiment { class Aging2Experiment; }
namespace gfe::experiment::details { class Aging2Worker; }
namespace gfe::experiment::details { class LatencyStatistics; }

namespace gfe::experiment::details {

class Aging2Master {
    friend class Aging2Worker;

    const Aging2Experiment& m_parameters;
    const bool m_is_directed; // is the graph directed?
    std::vector<Aging2Worker*> m_workers; // pool of workers
    std::atomic<int64_t> m_num_operations_performed = 0; // current number of operations performed so far
    std::atomic<int> m_last_progress_reported = 0; // the last progress of the experiment, reported by any of the worker threads. E.g. 1%, 2%, 3%, so on.

    // report how long it took to perform 1x, 2x, 3x, ... updates w.r.t. to the loaded graph.
    std::chrono::steady_clock::time_point m_time_start; // when the computation started
    uint64_t* m_reported_times = nullptr; // microsecs
    std::atomic<int> m_last_time_reported = 0;

    // latencies of each update
    uint64_t* m_latencies = nullptr; // nanosecs
    uint64_t m_latencies_num_insertions = 0; // total number of operations that are insertions
    uint64_t m_latencies_num_deletions = 0; // total number of operations that are deletions

    cuckoohash_map<uint64_t, bool> m_vertices_present; // current list of vertices present in the graph

    Aging2Result m_results; // final results of the experiment

    // Initialise the set of workers
    void init_workers();

    // Load & partition the edges to insert/remove in the available workers
    void load_edges();

    // Prepare the array to record the latency of all updates
    void prepare_latencies();

    // Execute the main part of the experiment, that is the insertions/deletions in the graph with the worker threads
    void do_run_experiment();

    // Remove the vertices that do not belong to the final graph
    void remove_vertices();

    // Save the current results in `m_results'
    void store_results();

    // Retrieve the current number of operations performed so far by the workers
    uint64_t num_operations_sofar() const;

    // Wait for the workers to complete, record the throughput in the meanwhile
    void wait_and_record();

    // print to stdout the number of vertices/edges expected and effectively stored in the final library. The two values should be equal.
    void log_num_vtx_edges();

    // Grab the vertex id of a random (final) edge
    void set_random_vertex_id(uint64_t* edges, uint64_t num_edges);

public:
    Aging2Master(const Aging2Experiment& parameters);

    // Destructor
    ~Aging2Master();

    // Execute the actual experiment
    Aging2Result execute();

    // Is the graph directed?
    bool is_directed() const { return m_is_directed; }

    // Total number of operations to perform (insertions/deletions)
    uint64_t num_operations_total() const;

    // Total number of edges expected in the final graph
    uint64_t num_edges_final_graph() const;

    // Access the configuration of this experiment
    const Aging2Experiment& parameters() const { return m_parameters; }
};

}
