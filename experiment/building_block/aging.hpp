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

#include "third-party/libcuckoo/cuckoohash_map.hh"

// forward declarations
namespace graph { class Edge; }
namespace graph { class WeightedEdge; }
namespace graph { class WeightedEdgeStream; }
namespace library { class UpdateInterface; }

namespace experiment {

class Aging {
    Aging(const Aging&) = delete;
    Aging& operator=(const Aging&) = delete;

    std::shared_ptr<library::UpdateInterface> m_interface; // the library where vertices and edges will be inserted
    std::shared_ptr<graph::WeightedEdgeStream> m_stream; // the final graph to insert

    cuckoohash_map<uint64_t, bool> m_vertices_present; // current list of vertices present
    enum EdgesPresentState { TEMPORARY, FINAL };
    cuckoohash_map<graph::Edge, EdgesPresentState> m_edges_present; // the current list of edges present in the map

    std::atomic<int64_t> m_num_operations_performed = 0; // current number of operations performed
    const uint64_t m_num_operations_total; // total number of operations to perform
    std::atomic<int64_t> m_startup_counter = 0; // keep track of the number of threads that still need to be initialised
    const int64_t m_num_threads; // the total number of concurrent threads to use

    double m_expansion_factor = 1.3; // set the maximum size, in terms of number of edges, that the underlying interface should contain

    uint64_t m_granularity = 1024; // the granularity


    // Insert the given edge
    void thread_insert(const graph::WeightedEdge& edge);

    // Execute the workload for the given thread
    void thread_main(int thread_id);

public:
    Aging(std::shared_ptr<library::UpdateInterface> interface);

    Aging(std::shared_ptr<library::UpdateInterface> interface, std::shared_ptr<graph::WeightedEdgeStream> stream);

    Aging(std::shared_ptr<library::UpdateInterface> interface, std::shared_ptr<graph::WeightedEdgeStream> stream, int64_t num_threads);

    // run the experiment
    std::chrono::microseconds execute();

    // Set the max expansion factor of the graph, in terms of number of active edges, for the underlying interface. It should be a value >= 1.
    // For instance with f = 1.3, then the graph can grow up to 30% than the size of the final graph
    void set_expansion_factor(double factor);

    // Set the granularity of a single round of operation inside a thread
    void set_operation_granularity(uint64_t granularity);
};

}

