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

#include <condition_variable>
#include <mutex>
#include <random>
#include <thread>

#include "common/circular_array.hpp"
#include "graph/edge.hpp"

// forward declarations
namespace experiment::details { class Aging2Master; }
namespace library { class UpdateInterface; }

namespace experiment::details {

class Aging2Worker {
    Aging2Worker(const Aging2Worker&) = delete;
    Aging2Worker& operator=(const Aging2Worker&) = delete;

    Aging2Master& m_master; // pointer to the master thread
    library::UpdateInterface* m_library; // the library being evaluated
    const int m_worker_id; // this id is passed to the interface #on_worker_init and #on_worker_destroy
    uint64_t m_partition_start; // the location where to insert the `real' edges that belong to the final graph
    uint64_t m_partition_end; // the last position (+1) of the partition of `real' edges to insert, from the final graph
    uint64_t m_partition_offset = 0; // offset from `m_partition_start' for the next to insert from the final graph
    common::CircularArray<graph::Edge> m_edges2remove; // edges that have been inserted but do not belong to the final graph
    std::mt19937_64 m_random { std::random_device{}() }; // pseudo-random generator
    std::uniform_real_distribution<double> m_uniform{ 0., 1. }; // uniform distribution in [0, 1]

    enum class Task { IDLE, START, STOP, EXECUTE_UPDATES, REMOVE_ARTIFICIAL_VERTICES };
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

    // execute the insert/delete operations for the graph in the background thread
    void main_execute_updates();

    // remote the artificial vertices, those that do not belong to the final graph, in the background thread
    void main_remove_artificial_vertices();

    // Set the task to execute asynchronously
    void set_task_async(Task t);

    // Insert a single vertex in the graph system/library, if it's not already present
    void graph_insert_vertex(uint64_t vertex_id);

    // Insert the given edge in the graph
    void graph_insert_edge(graph::WeightedEdge edge);

    // Remove the given edge from the graph
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

    // Request the thread to execute all updates
    void execute_updates();

    // Request to remove the vertices that do not belong to the final graph
    void remove_artificial_vertices();

    // Wait for the last operation issued to complete
    void wait();
};

} // namespace
