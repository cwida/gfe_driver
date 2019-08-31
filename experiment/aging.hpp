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

#include <atomic>
#include <chrono>
#include <cinttypes>
#include <memory>
#include <vector>

#include <gmpxx.h> // libgmp

#include "common/circular_array.hpp"
#include "details/aging_operation.hpp"
#include "details/aging_partition.hpp"
#include "graph/edge.hpp"
#include "third-party/libcuckoo/cuckoohash_map.hh"

// forward declarations
namespace graph { class WeightedEdgeStream; }
namespace experiment::details { class AgingThread; }
namespace experiment::details { class AsyncBatch; }
namespace library { class UpdateInterface; }

namespace experiment {

class Aging {
    friend class details::AgingThread;

    Aging(const Aging&) = delete;
    Aging& operator=(const Aging&) = delete;

    std::shared_ptr<library::UpdateInterface> m_interface; // the library where vertices and edges will be inserted
    std::shared_ptr<graph::WeightedEdgeStream> m_stream; // the final graph to insert

    cuckoohash_map<uint64_t, bool> m_vertices_present; // current list of vertices present
    const cuckoohash_map<uint64_t, bool> m_vertices_final; // list of vertices in the final graph

    std::atomic<int64_t> m_num_operations_performed = 0; // current number of operations performed
    const uint64_t m_num_operations_total; // total number of operations to perform
    const int64_t m_num_threads; // the total number of concurrent threads to use
    const uint64_t m_num_edges; // the total number of edges in the final graph
    const uint64_t m_num_vertices; // the total number of vertices in the final graph (a bit expensive to compute, but ok, c'est la vie)
    uint64_t m_max_vertex_id_artificial; // the maximum vertex id that can be generated
    const uint64_t m_max_vertex_id_final; // the maximum vertex id in the final graph
    const bool m_is_directed; // whether the graph is directed
    const double m_max_weight; // the max weight for an edge in graph, incl.
    uint64_t* m_vertices2remove = nullptr; // at the end of the experiment, it's an array with the whole amount of vertices in excess that need to be removed from the final graph

    double m_ef_edges = 1; // multiplying constant on the max number of edges, stored at a given time, w.r.t. to the loaded
    double m_ef_vertices = 1; // multiplying constant on the max number of vertices that can be generated (const = 1, "almost" no new vertices)
    int64_t m_granularity = 1024; // the granularity of each burst of insertions/deletions, as executed by a worker thread
    uint64_t m_batch_size = 0; // send the updates in batches
    std::chrono::milliseconds m_build_frequency {0}; // the frequency to create a new delta/snapshot, that is invoking the method #build()

    // final data
    uint64_t m_num_artificial_vertices = 0; // the total number of artificial vertices (not present in the loaded graph)  inserted in the updates
    uint64_t m_completion_time = 0; // the amount of time to complete all updates, in microsecs
    uint64_t m_num_vertices_final_graph = 0; // the number of vertices in the final graph, after all updates have been performed
    uint64_t m_num_edges_final_graph = 0; // the number of edges in the final graph, after all updates have been performed
    uint64_t m_num_build_invocations = 0; // total number of invocations to the method #build

    // whether to report the current progress
    bool m_report_progress = false;
    std::atomic<int> m_last_progress_reported = 0;

    // report how long it took to perform 1x, 2x, 3x, ... updates w.r.t. to the loaded graph.
    std::chrono::steady_clock::time_point m_time_start; // when the computation started
    uint64_t* m_reported_times = nullptr; // microsecs
    std::atomic<int> m_last_time_reported = 0;

    // Insert a single vertex in the graph system/library, if it's not already present. Invoked concurrently by the worker threads
    void insert_vertex(uint64_t vertex_id);

    // Execute the given action on all workers
    void all_workers_execute(const std::vector<std::unique_ptr<details::AgingThread>>& workers, details::AgingOperation operation);

    // print to stdout the number of vertices/edges expected and effectively stored in the final library. The two values should be equal.
    void log_num_vtx_edges();

    std::vector<std::vector<details::AgingPartition>> compute_partitions_undirected() const;
    std::vector<std::vector<details::AgingPartition>> compute_partitions_directed() const;
    void compute_partitions_add_missing(std::vector<std::vector<details::AgingPartition>>* out_partitions) const;


public:
    /**
     * Constructor
     * @param interface the library to evaluate
     * @param stream reader for the final graph to load
     * @param mult_num_operations coefficient for the total number of operations to perform, the final amount is given by mult_num_operations * stream->num_edges()
     * @param num_threads number of threads to use
     */
    Aging(std::shared_ptr<library::UpdateInterface> interface, std::shared_ptr<graph::WeightedEdgeStream> stream, double mult_num_operations, int64_t num_threads);


    Aging(std::shared_ptr<library::UpdateInterface> interface, std::shared_ptr<graph::WeightedEdgeStream> stream, uint64_t num_operations, int64_t num_threads, bool is_directed, double max_weight);

    // destructor
    ~Aging();

    // run the experiment
    std::chrono::microseconds execute();

    // check whether the graph is directed
    bool is_directed() const;

    // Set the max expansion factor of the graph, related to the number of active edges, for the underlying interface. It should be a value >= 1.
    // For instance with f = 1.3, then the graph can grow up to 30% than the size of the final graph
    void set_expansion_factor_edges(double factor);

    // Set the max expansion factor of the graph, related to the number of vertices. It should be a value >= 1.
    // For instance with f = 1.3, then the experiment should create 30% more artificial vertices than those expected in the final graph
    void set_expansion_factor_vertices(double factor);

    // Set the granularity of a single round of operation inside a thread
    void set_operation_granularity(uint64_t granularity);

    // Request to send edge updates in batches of the given size. Vertex insertions will continue to be sent one at the time
    void set_batch_size(uint64_t size);

    // Whether to print to stdout the current progress
    void set_report_progress(bool value);

    // Set how frequently create a new snapshot/delta in the library (0 = do not create new snapshots)
    void set_build_frequency(std::chrono::milliseconds millisecs);

    // Store the results into the database
    void save();
};

}

