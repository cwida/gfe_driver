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
#include <utility>
#include "reader.hpp"

namespace graph { class WeightedEdge; } // forward decl.

namespace reader {

/**
 * Loader for the METIS v5.1 file format
 */
class MetisReader : public Reader {
    std::fstream m_handle; // internal file handle
    std::string m_current_line; // the current line being processed
    char* m_current_line_read_value_ptr { nullptr }; // the last position parsed by #read_value();
    uint64_t m_num_vertices { 0 }; // the number of vertices present in the graph, according to the file header
//    uint64_t m_num_edges { 0 }; // the number of edges present in the graph, according to the file header (an undirected edge is counted only once)
    int m_num_vertex_weights { 0 }; // the number of vertex weights + extra field `size', according to the file header. These fields are skipped when parsing the file
    bool m_has_edge_weight { false }; // whether each edge contains a weight
    bool m_is_edge_weight_constant { false }; // whether to assume that all edge weights are not provided and implicitly set to 1
    std::mt19937_64 m_random_generator; // random generator for non weighted graphs, it will assign a random weight, used when both m_has_edge_weight == false && m_is_edge_weight_constant == false

    // current line being parsed
    uint64_t m_lineno = 0;

    // The last edge parsed from the input file
    uint64_t m_edge_vertex1 { 0 }, m_edge_vertex2 { 0 };
    uint32_t m_edge_weight { 0 };

    /**
     * Check whether the current line is a comment (the line starts with %)
     */
    bool is_comment() const;

    /**
     * Fetch the next line which is not a comment
     */
    bool fetch_next_line();

    /**
     * Parse the next uint64_t value from the current line
     */
    std::pair<bool, uint64_t> read_value();

public:
    /**
     * Read the edge list from the given file
     */
    MetisReader(const std::string& path);

    /**
     * Destructor
     */
    ~MetisReader();

    /**
     * Interface, report one edge at the time
     */
    bool read(graph::WeightedEdge&) override;

    /**
     * Interface, whether the input graph is directed
     */
    bool is_directed() const override;
};


} // namespace reader





