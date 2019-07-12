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

#include "reader.hpp"

#include <unordered_map>

namespace graph { class WeightedEdge; } // forward decl.

namespace reader {

/**
 * Parser to read the edge list from the LDBC benchmarks
 */
class GraphalyticsReader : public Reader {
    std::unordered_map<std::string, std::string> m_properties; // property file
    bool m_directed = true; // whether the graph being processed is directed or not
    bool m_is_weighted = true; // whether the graph being processed contains weights or not
    void* m_handle_edge_file { nullptr }; // fstream, handle to parse the edge-file
    void* m_handle_vertex_file { nullptr }; // fstream, handle to parse the vertex-file
    uint64_t m_last_source {0}; uint64_t m_last_destination {0}; double m_last_weight {0.0}; // the last edge being parsed
    bool m_last_reported = true; // whether we have reported the last edge with source/dest vertices swapped in an undirected graph
    bool m_emit_directed_edges = false; // if the graph is undirected, report the same edge twice as src -> dest and dest -> src

    // Close the internal handles to parse the edge-file/vertex-file
    void close();

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
     * Init the reader with the path to the graph property files (*.properties)
     */
    GraphalyticsReader(const std::string& path_properties);

    /**
     * Destructor
     */
    ~GraphalyticsReader();

    /**
     * Interface, report one edge at the time
     */
    bool read(graph::WeightedEdge&) override;

    /**
     * Iterator, read one edge at the time from the graph
     */
    bool read_edge(graph::WeightedEdge& out_edge);

    /**
     * Iterator, read one vertex at the time from the graph
     */
    bool read_vertex(uint64_t& out_vertex);

    /**
     * Reset the position of the iterators read/read_edge/read_vertex at the start of the file
     */
    void reset();

    // Retrieve the given property from the map, or the empty string if the property is not present
    std::string get_property(const std::string& key) const;

    /**
     * Path to the vertex file
     */
    std::string get_path_vertex_list() const;

    /**
     * Path to the edge file
     */
    std::string get_path_edge_list() const;

    /**
     * Check whether the graph is directed
     */
    bool is_directed() const override;

    /**
     * Check whether the graph is weighted
     */
    bool is_weighted() const;

    /**
     * Shall we report the same edge twice in an undirected graph, once as src -> dest and once as dest -> src
     */
    void set_emit_directed_edges(bool value);
};


} // namespace


