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
#include <condition_variable>
#include <mutex>
#include <random>
#include <thread>

#include "common/circular_array.hpp"
#include "graph/edge.hpp"

// forward declarations
namespace gfe::experiment::details { class Aging2Master; }
namespace gfe::library { class UpdateInterface; }

namespace gfe::experiment::details {

class Aging2Worker {
    Aging2Worker(const Aging2Worker&) = delete;
    Aging2Worker& operator=(const Aging2Worker&) = delete;

    Aging2Master& m_master; // pointer to the master thread
    library::UpdateInterface* m_library; // the library being evaluated
    const int m_worker_id; // this id is passed to the interface #on_worker_init and #on_worker_destroy
    common::CircularArray<std::vector<gfe::graph::WeightedEdge>*> m_updates; // the updates to perform
    std::mt19937_64 m_random { std::random_device{}() }; // pseudo-random generator
    std::uniform_real_distribution<double> m_uniform{ 0., 1. }; // uniform distribution in [0, 1]
    uint64_t* m_latency_insertions {nullptr};
    uint64_t m_num_edge_insertions {0}; // counter, total number of edge insertions to perform, as contained in the array m_updates
    uint64_t* m_latency_deletions {nullptr};
    uint64_t m_num_edge_deletions {0}; // counter, total number of edge deletions to perform, as contained in the array m_updates
    std::atomic<uint64_t> m_num_operations = 0; // counter, total number of operations performed so far

    enum class TaskOp { IDLE, START, STOP, LOAD_EDGES, EXECUTE_UPDATES, REMOVE_VERTICES, SET_ARRAY_LATENCIES };
    struct Task { TaskOp m_type; uint64_t* m_payload; uint64_t m_payload_sz; };
    Task m_task; // current task being executed

    std::thread m_thread; // the thread associated to the background thread/service
    std::mutex m_mutex;
    std::condition_variable m_condvar;

    // start the service, i.e. the execution of the background thread
    void start();

    // stop the service, i.e. the execution of the background thread
    void stop();

    // the controller for the background thread
    void main_thread();

    // load a batch of edges in the background thread
    void main_load_edges(uint64_t* edges, uint64_t num_edges);

    // execute the insert/delete operations for the graph in the background thread
    void main_execute_updates();

    // remote the artificial vertices, those that do not belong to the final graph, in the background thread
    void main_remove_vertices(uint64_t* vertices, uint64_t num_vertices);

    // Set the task to execute asynchronously
    void set_task_async(TaskOp task_type, uint64_t* payload = nullptr, uint64_t payload_sz = 0);

    // Execute a batch of updates
    void graph_execute_batch_updates(graph::WeightedEdge* __restrict updates, uint64_t num_updates);

    template<bool with_latency>
    void graph_execute_batch_updates0(graph::WeightedEdge* __restrict updates, uint64_t num_updates);

    // Insert a single vertex in the graph system/library, if it's not already present
    void graph_insert_vertex(uint64_t vertex_id);

    // Insert the given edge in the graph
    template<bool with_latency>
    void graph_insert_edge(graph::WeightedEdge edge);

    // Remove the given edge from the graph
    template<bool with_latency>
    void graph_remove_edge(graph::Edge edge, bool force = true);

    // Remove the temporary edge at the head of the queue m_edges2remove
    void graph_remove_temporary_edge();

    // The size of each burst of insertions/deletions
    uint64_t granularity() const;

    // How many edges do we still of the final graph do we still need to insert
    int64_t missing_edges_final() const;

public:
    Aging2Worker(Aging2Master& master, int worker_id);

    // Destructor. It implicitly stops the background thread.
    ~Aging2Worker();

    // Load a batch of edges
    void load_edges(uint64_t* edges, uint64_t num_edges);

    // Set the latency arrays for insertions and deletions
    void set_latencies(uint64_t* array_insertions, uint64_t* array_deletions);

    // Request the thread to execute all updates
    void execute_updates();

    // Request to remove the vertices that do not belong to the final graph
    void remove_vertices(uint64_t* vertices, uint64_t num_vertices);

    // Wait for the last operation issued to complete
    void wait();

    // Wait for the last operation issued to complete or up to the given time point. Return true if the last operation performed was completed.
    bool wait(const std::chrono::time_point<std::chrono::steady_clock>& tp);

    // Number of edge insertions scheduled to perform
    uint64_t num_insertions() const;

    // Number of edge deletions scheduled to perform
    uint64_t num_deletions() const;

    // Total number of operations performed so far
    uint64_t num_operations() const;
};

} // namespace
