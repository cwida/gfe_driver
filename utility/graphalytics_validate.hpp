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

#include <string>

#include "common/error.hpp"

namespace utility {

// In case the result does not match the reference result, the should should throw a GraphalyticsValidateError
DEFINE_EXCEPTION(GraphalyticsValidateError);

/**
 * Validate the result of an algorithm from the Graphalytics interface with its reference/expected output.
 */
class GraphalyticsValidate {
protected:
    // The two files should be identical
    static void exact_match(const std::string& result, const std::string& expected);

public:

    /**
     * Validate the output of the BFS (result) with the given reference (expected)
     * @throw GraphalyticsValidateError in case of mismatch
     */
    static void bfs(const std::string& result, const std::string& expected);

};
} // namespace


