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

#include <cinttypes>
#include <string>
#include <unordered_map>

#include "common/error.hpp"

namespace gfe::utility {

// In case the result does not match the reference result, the should should throw a GraphalyticsValidateError
DEFINE_EXCEPTION(GraphalyticsValidateError);

/**
 * Validate the result of an algorithm from the Graphalytics interface with its reference/expected output.
 */
class GraphalyticsValidate {
public:
    using vertex_map_t = std::unordered_map<uint64_t, uint64_t>; // remap the vertices from the `expected file' into the `result file'

protected:
    // The two files should be identical
    static void exact_match(const std::string& result, const std::string& expected, uint64_t max_num_errors, const vertex_map_t* vtx_map, bool vtx_relabel_values);

    // The two files contain values of type double, the difference between r (result) and e (expected) should be small enough: |r - e| / e < epsilon, with epsilon parameter
    static void epsilon_match(const std::string& result, const std::string& expected, double epsilon, uint64_t max_num_errors, const vertex_map_t* vtx_map);

    // There must exist a bijective function f to map the values val from the two files, i.e. f[ val(v1) ] = val(v2)
    static void equivalence_match(const std::string& result, const std::string& expected, uint64_t max_num_errors, const vertex_map_t* vertex_map);

public:

    /**
     * Validate the output of the BFS (result) with the given reference (expected)
     * @throw GraphalyticsValidateError in case of mismatch
     */
    static void bfs(const std::string& result, const std::string& expected, uint64_t max_num_errors = 1, const vertex_map_t* vertex_map = nullptr);

    /**
     * Validate the output of the PageRank algorithm (result) with the given reference (expected)
     * @throw GraphalyticsValidateError in case of mismatch
     */
    static void pagerank(const std::string& result, const std::string& expected, uint64_t max_num_errors = 1, const vertex_map_t* vertex_map = nullptr);

    /**
     * Validate the output of the WCC algorithm (result) with the given reference (expected)
     * @throw GraphalyticsValidateError in case of mismatch
     */
    static void wcc(const std::string& result, const std::string& expected, uint64_t max_num_errors = 1, const vertex_map_t* vertex_map = nullptr);

    /**
     * Validate the output of the LCC algorithm (result) with the given reference (expected)
     * @throw GraphalyticsValidateError in case of mismatch
     */
    static void lcc(const std::string& result, const std::string& expected, uint64_t max_num_errors = 1, const vertex_map_t* vertex_map = nullptr);

    /**
     * Validate the output of the CDLP algorithm (result) with the given reference (expected)
     * @throw GraphalyticsValidateError in case of mismatch
     */
    static void cdlp(const std::string& result, const std::string& expected, uint64_t max_num_errors = 1, const vertex_map_t* vertex_map = nullptr);

    /**
     * Validate the output of the SSSP algorithm (result) with the given reference (expected)
     * @throw GraphalyticsValidateError in case of mismatch
     */
    static void sssp(const std::string& result, const std::string& expected, uint64_t max_num_errors = 1, const vertex_map_t* vertex_map = nullptr);

};

} // namespace
