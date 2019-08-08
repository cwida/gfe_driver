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
#include <ostream>
#include <string>
#include <vector>

namespace library { class GraphalyticsInterface; } // forward decl.

namespace experiment {

/**
 * The properties of the algorithm to run in the Graphalytics suite
 */
struct GraphalyticsAlgorithms {
    struct {
        bool m_enabled = false;
        uint64_t m_source_vertex = 0;
    } bfs;
    struct {
        bool m_enabled = false;
        uint64_t m_max_iterations = 0;
    } cdlp;
    struct {
        bool m_enabled = false;
    } lcc;
    struct {
        bool m_enabled = false;
        double m_damping_factor = 0.85; // default
        double m_num_iterations = 0;
    } pagerank;
    struct {
        bool m_enabled = false;
        uint64_t m_source_vertex = 0;
    } sssp;
    struct {
        bool m_enabled = false;
    } wcc;


    /**
     * Load the properties of each algorithms from the given files
     */
    GraphalyticsAlgorithms(const std::string& path);
};

std::ostream& operator<<(std::ostream& out, const GraphalyticsAlgorithms& props); // debug only


/**
 * Execute one by one the algorithms of the graphalytics suite, up to N times
 */
class GraphalyticsSequential{
    std::shared_ptr<library::GraphalyticsInterface> m_interface; // the library to evaluate
    const uint64_t m_num_repetitions; // number of times to repeat the execution of each algorithm
    const GraphalyticsAlgorithms m_properties; // the properties of the graphalytics algorithms

    // the completion times for each execution
    std::vector<int64_t> m_exec_bfs;
    std::vector<int64_t> m_exec_cdlp;
    std::vector<int64_t> m_exec_lcc;
    std::vector<int64_t> m_exec_pagerank;
    std::vector<int64_t> m_exec_sssp;
    std::vector<int64_t> m_exec_wcc;


public:
    GraphalyticsSequential(std::shared_ptr<library::GraphalyticsInterface> interface, uint64_t num_repetitions, const GraphalyticsAlgorithms& properties);

    /**
     * Execute the experiment
     */
    std::chrono::microseconds execute();

    /**
     * Report the execution results
     * @param save_in_db: if true store the results in the database
     */
    void report(bool save_in_db);
};


} // namespace exp
