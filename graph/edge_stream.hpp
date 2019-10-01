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

#include "edge.hpp"

#include <memory>
#include <string>
#include <vector>
#include "third-party/libcuckoo/cuckoohash_map.hh"


namespace graph {

class CByteArray; // forward decl.
class VertexList; // forward decl.

// Load in memory a stream of edges with
class WeightedEdgeStream {
    WeightedEdgeStream(const WeightedEdgeStream&) = delete;
    WeightedEdgeStream& operator=(const WeightedEdgeStream&) = delete;

    // edges are stores in three separate arrays:
    CByteArray* m_sources { nullptr }; // one for the sources
    CByteArray* m_destinations { nullptr }; // one for the destinations
    std::vector<double> m_weights; // and one for the weights

    // total number of edges
    uint64_t m_num_edges { 0 };

    // keep track of the maximum values for the vertices and the weights
    uint64_t m_max_vertex_id { 0 };
    double m_max_weight { 0 };

public:
    /**
     * Load the list of edges from the given file
     */
    WeightedEdgeStream(const std::string& path);

    /**
     * Load the list of edges from the given vector
     */
    WeightedEdgeStream(const std::vector<WeightedEdge>& vector);

    /**
     * Destructor
     */
    ~WeightedEdgeStream();

    // Perform a random permutation of the edge list. Use the default random seed.
    void permute();

    // Perform a random permutation of the edge list. The permutation is based on the given random seed.
    void permute(uint64_t seed);

    // Retrieve the edge at the given position, in [0, num_edges() )
    WeightedEdge get(uint64_t index) const;
    WeightedEdge operator[](uint64_t index) const { return get(index); } // alias

    // Retrieve the total number of edges in the list
    uint64_t num_edges() const { return m_num_edges; }

    // Get the list of vertices present
    std::unique_ptr<VertexList> vertex_list() const;

    // Get the set of vertices present. For each vertex, return the number of outgoing edges attached.
    std::unique_ptr<cuckoohash_map<uint64_t, uint64_t>> vertex_table() const;

    // Get the max vertex id present
    uint64_t max_vertex_id() const { return m_max_vertex_id; };

    // Get the max weight for an edge in the graph
    double max_weight() const { return m_max_weight; }
};

} // namespace graph

