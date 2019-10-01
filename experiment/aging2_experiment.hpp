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
#include <cstdint>
#include <memory>

#include "aging2_result.hpp"

// forward declarations
namespace graph { class WeightedEdgeStream; }
namespace experiment { class Aging2Experiment; }
namespace experiment::details { class Aging2Master; }
namespace experiment::details { class Aging2Worker; }
namespace library { class UpdateInterface; }

namespace experiment {

/**
 * Builder/factory class to create & execute instances of the Aging experiment.
 *
 * This class is not thread-safe.
 */
class Aging2Experiment {
    Aging2Experiment(const Aging2Experiment&) = delete;
    Aging2Experiment& operator=(const Aging2Experiment&) = delete;
    friend class Aging2Result;
    friend class details::Aging2Master;
    friend class details::Aging2Worker;

    std::shared_ptr<library::UpdateInterface> m_library; // the library to evaluate
    std::shared_ptr<graph::WeightedEdgeStream> m_stream; // the graph to simulate. After the experiment, the library will contain the same nodes/edges of the given graph.
    uint64_t m_num_threads = 1; // set the number of threads to use
    uint64_t m_worker_granularity = 1024; // the granularity of a task for a worker, that is the number of contiguous operations (inserts/deletes) performed inside the threads between each invocation to the scheduler.
    double m_mult_ops = 10.0; // the number of operations to perform, w.r.t. to the size of the provided graph. It must be a value >= 1.
    double m_ef_edges = 1; // multiplying constant on the max number of edges, stored at a given time, w.r.t. to the loaded
    double m_ef_vertices = 1; // multiplying constant on the max number of vertices that can be generated (const = 1, "almost" no new vertices)
    double m_max_weight = 1024.0; // set the max weight for the edges to create
    std::chrono::milliseconds m_build_frequency {0}; // the frequency to create a new delta/snapshot, that is invoking the method #build()
    bool m_report_progress = false; // whether to report the current progress

public:
    // Instantiate the factory class
    Aging2Experiment();

    // Set the library to evaluate
    void set_library(std::shared_ptr<library::UpdateInterface> library);

    // Set the graph to produce
    void set_graph(std::shared_ptr<graph::WeightedEdgeStream> graph);

    // Set the number of operations (inserts/deletes) to perform, w.r.t. to the size of the given graph
    void set_coeff_operations(double coeff_operations);

    // Set the max weight for the edges created
    void set_max_weight(double value);

    // Set the number of client threads to use in the experiment, that is, the parallelism degree
    void set_parallelism_degree(uint64_t num_threads);

    // Set the max expansion factor of the graph, related to the number of active edges, for the underlying interface. It should be a value >= 1.
    // For instance with f = 1.3, then the graph can grow up to 30% than the size of the final graph
    void set_expansion_factor_edges(double factor);

    // Set the max expansion factor of the graph, related to the number of vertices. It should be a value >= 1.
    // For instance with f = 1.3, then the experiment should create 30% more artificial vertices than those expected in the final graph
    void set_expansion_factor_vertices(double factor);

    // Set how frequently create a new snapshot/delta in the library (0 = do not create new snapshots)
    void set_build_frequency(std::chrono::milliseconds millisecs);

    // Whether to print to stdout the current progress of the experiment
    void set_report_progress(bool value);

    // [Internal parameter]
    // Set the granularity of a task for a worker thread. This is the number of contiguos operations (inserts/deletes) done
    // by each worker thread between each invocation to the scheduler.
    void set_worker_granularity(uint64_t value);

    // Execute the experiment with the given configuration
    // @param reset_graph if true, release the contained graph before running the experiment, to save some memory
    Aging2Result execute();
};

} // namespace
