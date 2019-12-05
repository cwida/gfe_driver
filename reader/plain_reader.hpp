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

#include <fstream>
#include <random>
#include "reader.hpp"

namespace gfe::graph { class WeightedEdge; } // forward decl.

namespace gfe::reader {

/**
 * Specific implementation to load the edges from a plain representation (.el or .wel). The file
 * should contain an edge per line:
 * - for non weighted graphs, a line in the format: src dst
 * - for weighted graphs, a line in the format: src dst weight
 * Any line starting with a sharp sign (#) is considered a comment and it is ignored
 */
class PlainReader : public Reader {
    std::fstream m_handle; // internal file handle
    const bool m_is_weighted; // whether we are reading a weighted graph
    const uint32_t m_max_weight; // for non weighted graphs
    std::mt19937_64 m_random_generator; // used for non weighted graphs to generate a random weight in [1, m_max_weight]

    /**
     * Check whether the given line is a comment or empty, that is, it starts with a sharp symbol # or contains no symbols
     */
    static bool ignore_line(const char* line);
    static bool ignore_line(const std::string& line);

    /**
     * Check whether the current marker points to a number
     */
    static bool is_number(const char* marker);

public:
    /**
     * Read the edge list from the given file
     * @param path the source of the file
     * @param is_weighted whether the input file is weighted (file.wel) or not (file.el)
     */
    PlainReader(const std::string& path, bool is_weighted);

    /**
     * Read the edge list from the given file
     * @param path the source of the file
     * @param is_weighted whether the input file is weighted
     * @param max_weight if the graph is non weighted, then it's the maximum weight that can be assigned, at random, to each parsed edge
     */
    PlainReader(const std::string& path, bool is_weighted, uint32_t max_weight);

    /**
     * Destructor
     */
    ~PlainReader();

    // interface
    bool read(graph::WeightedEdge&) override;
    bool is_directed() const override;
};

} // namespace
