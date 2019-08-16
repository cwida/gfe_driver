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
#include <condition_variable>
#include <cinttypes>
#include <future>
#include <memory>
#include <mutex>
#include <random>
#include <unordered_map>
#include <vector>

#include <gmpxx.h> // libgmp

#include "common/circular_array.hpp"
#include "graph/edge.hpp"
#include "library/interface.hpp"
#include "third-party/libcuckoo/cuckoohash_map.hh"

// forward declarations
namespace graph { class Edge; }
namespace graph { class WeightedEdge; }
namespace graph { class WeightedEdgeStream; }
namespace experiment::details { class AsyncBatch; }

namespace experiment {

class Aging {
    Aging(const Aging&) = delete;
    Aging& operator=(const Aging&) = delete;

    std::shared_ptr<library::UpdateInterface> m_interface; // the library where vertices and edges will be inserted
    std::shared_ptr<graph::WeightedEdgeStream> m_stream; // the final graph to insert

    cuckoohash_map<uint64_t, bool> m_vertices_present; // current list of vertices present
    const cuckoohash_map<uint64_t, bool> m_vertices_final; // list of vertices in the final graph

    std::atomic<int64_t> m_num_operations_performed = 0; // current number of operations performed
    const uint64_t m_num_operations_total; // total number of operations to perform
    const int64_t m_num_threads; // the total number of concurrent threads to use
    const int64_t m_num_edges; // the total number of edges in the final graph
    const uint64_t m_max_vertex_id; // the maximum vertex id that can be generated
    const bool m_is_directed; // whether the graph is directed
    const double m_max_weight; // the max weight for an edge in graph, incl.

    double m_expansion_factor = 1.3; // set the maximum size, in terms of number of edges, that the underlying interface should contain
    int64_t m_granularity = 1024; // the granularity of each burst of insertions/deletions, as execute by a worker thread
    uint64_t m_batch_size = 0; // send the updates in batches

    // A single partition handled by a worker thread
    // Because these are vertex_ids, these can be stored in 64 bits
    struct Partition { uint64_t m_start, m_length; Partition(uint64_t s, uint64_t l) : m_start(s), m_length(l) { }; };

    // The single operations that can be performed by the worker threads
    enum class AgingOperation { NONE, START, STOP, COMPUTE_FINAL_EDGES, EXECUTE_EXPERIMENT, INTERNAL_CLEANUP, REMOVE_VERTICES };

    uint64_t m_completion_time = 0; // the amount of time to complete all updates, in microsecs

    uint64_t* m_vertices2remove = nullptr; // at the end of the experiment, it's an array with the whole amount of vertices in excess that need to be removed from the final graph

    // A single worker thread in the aging experiment
    class AgingThread {
        Aging* m_instance; // aging instance
        library::UpdateInterface* m_interface; // the graph library we are operating on
        const int m_worker_id; // the id of this thread
        const bool m_is_undirected; // is the graph undirected?

        std::vector<Partition> m_partitions;
        uint64_t m_num_src_vertices_in_partitions {0};

        std::vector<graph::WeightedEdge> m_edges; // primary edges to insert in the final graph
        uint64_t m_final_edges_current_position { 0 }; // index to keep track up to where which edges have been inserted

        std::unordered_map<graph::Edge, bool> m_edges_already_inserted; // keep track of which edges has already been inserted
        common::CircularArray<graph::Edge> m_edges2remove; // edges that have been inserted but do not belong to the final graph

        experiment::details::AsyncBatch* m_batch; // send the edge updates in batches

        // synchronisation with the master thread
        AgingOperation m_current_operation { AgingOperation::NONE }; // the current operation the thread is performing
        std::promise<void> m_callback;
        std::mutex m_mutex_op;
        std::condition_variable m_condvar_op;

        std::mt19937_64 m_random { std::random_device{}() };
        std::uniform_real_distribution<double> m_uniform{ 0., 1. };

        // arbitrary precision arithmetic
        mpz_class* m_num_edges_in_partition { nullptr } ; // array of size m_partitions.size(), containing the max number of edges per partition

        Partition m_interval_vertices2remove { 0, 0 };

        // Transform a relative source id (randomly generated) into an absolute edge id, according to the handled partitions
        uint64_t src_rel2abs(uint64_t relative_vertex_id) const;

        // How many edges do we still of the final graph do we still need to insert
        int64_t missing_edges_final() const;

        // Insert the given edge
        void insert_edge(graph::WeightedEdge edge);

        // Remove the given edge
        void remove_edge(graph::Edge edge);

        // Remove the temporary edge at the head of the queue m_edges2remove
        void remove_temporary_edge();

        // Total number of edges that can be generated by this worker, assuming the graph is undirected
        mpz_class get_num_edges_in_my_partitions() const;

        // Transform the given edge id in the vertex ids src/dst
        void edge_id_2_vertices_id(const mpz_class& edge_id, uint64_t* out_src_id, uint64_t* out_dst_id);

        // Check whether the given vertex id (src) belongs to the set of partitions to handle
        bool vertex_belongs(uint64_t vertex_id) const;

        // Generate a random number in [0, 1]
        double random01() noexcept;

        // Is the graph undirected?
        bool is_undirected() const { return m_is_undirected; }

        // Flush the current batch of updates & reset the variable batch_pos to 0
        void batch_flush();

        // Logic to run the experiment
        void main_experiment();

    public:
        // Initialise the worker with the list of partitions in the adjacency matrix to handle
        AgingThread(Aging* instance, const std::vector<Partition>& partitions, int worker_id);

        // Destructor
        ~AgingThread();

        // Execute a single operation in the thread
        std::future<void> execute(AgingOperation operation);

        // Controller thread, sync with the caller
        void main_thread();

        // The last pass for this class is REMOVE_VERTICES, which requires to remove all vertices from the interface/library as set in the
        // array m_instance->m_vertices2remove in the given interval [start, start +length).
        void set_partition_vertices_to_remove(uint64_t start, uint64_t length);

    };
    friend class AgingThread;

    // Insert a single vertex in the graph system/library. Invoked concurrently by the worker threads
//    void insert_edge(graph::WeightedEdge edge);

    // Insert a single vertex in the graph system/library, if it's not already present. Invoked concurrently by the worker threads
    void insert_vertex(uint64_t vertex_id);

    // Execute the given action on all workers
    void all_workers_execute(const std::vector<std::unique_ptr<AgingThread>>& workers, AgingOperation operation);

    std::vector<std::vector<Partition>> compute_partitions_undirected() const;
    std::vector<std::vector<Partition>> compute_partitions_directed() const;

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

    // run the experiment
    std::chrono::microseconds execute();

    // check whether the graph is directed
    bool is_directed() const;

    // Set the max expansion factor of the graph, in terms of number of active edges, for the underlying interface. It should be a value >= 1.
    // For instance with f = 1.3, then the graph can grow up to 30% than the size of the final graph
    void set_expansion_factor(double factor);

    // Set the granularity of a single round of operation inside a thread
    void set_operation_granularity(uint64_t granularity);

    // Request to send edge updates in batches of the given size. Vertex insertions will continue to be sent one at the time
    void set_batch_size(uint64_t size);

    // Store the results into the database
    void save();
};

}

