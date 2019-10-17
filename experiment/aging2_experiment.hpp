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
    friend class Aging2Result;
    friend class details::Aging2Master;
    friend class details::Aging2Worker;

    std::shared_ptr<library::UpdateInterface> m_library; // the library to evaluate
    std::string m_path_log; // the path to the log file [graphlog] with the sequence of updates to perform
    uint64_t m_num_threads = 1; // set the number of threads to use
    uint64_t m_worker_granularity = 1024; // the granularity of a task for a worker, that is the number of contiguous operations (inserts/deletes) performed inside the threads between each invocation to the scheduler.
    double m_max_weight = 1024.0; // set the max weight for the edges to create
    std::chrono::milliseconds m_build_frequency {0}; // the frequency to create a new delta/snapshot, that is invoking the method #build()
    bool m_report_progress = false; // whether to report the current progress

public:
    // Instantiate the factory class
    Aging2Experiment();

    // Set the library to evaluate
    void set_library(std::shared_ptr<library::UpdateInterface> library);

    // Set the path to the log file with all updates
    void set_log(const std::string& path_log);

    // Set the max weight for the edges created
    void set_max_weight(double value);

    // Set the number of client threads to use in the experiment, that is, the parallelism degree
    void set_parallelism_degree(uint64_t num_threads);

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
