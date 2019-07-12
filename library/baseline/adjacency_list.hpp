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

#include "common/error.hpp"
#include "library/interface.hpp"

#include <cinttypes>
#include <mutex>
#include <unordered_map>
#include <vector>

namespace library {

// Generic exception thrown by this class
DEFINE_EXCEPTION(AdjacencyListError);

/**
 * Sequential and base implementation of the interface, for testing purposes.
 * The class is thread-safe, but all exposed operations are serialised and sequential.
 */
class AdjacencyList : public virtual UpdateInterface, public virtual LoaderInterface {
    using EdgeList = std::vector</* edge */ std::pair< /* destination */ uint64_t,  /* weight */ double>>;
    using EdgePair = std::pair< /* outgoing edges */ EdgeList, /* incoming edges */ EdgeList >;
    using NodeList = std::unordered_map</* src vertex */ uint64_t, /* outgoing & incoming edges */ EdgePair>;
    NodeList m_adjacency_list;
    uint64_t m_num_edges = 0; // number of directed edges
    const bool m_is_directed; // whether the graph is directed or not

    using mutex_t = std::recursive_mutex;
    mutable mutex_t m_mutex; // read-write mutex

public:
    /**
     * Initialise the graph instance
     */
    AdjacencyList(bool is_directed);

    /**
     * Destructor
     */
    ~AdjacencyList();

    /**
     * Get the number of edges contained in the graph
     */
    virtual uint64_t num_edges() const;

    /**
     * Get the number of nodes stored in the graph
     */
    virtual uint64_t num_vertices() const;

    /**
     * Is the graph directed?
     */
    virtual bool is_directed() const;

    /**
     * Returns true if the given vertex is present, false otherwise
     */
    virtual bool has_vertex(uint64_t vertex_id) const;

    /**
     * Retrieve the weight associated to the given edge, or -1 if the given edge does not exist
     */
    virtual double get_weight(uint64_t source, uint64_t destination) const;

    /**
     * Dump the content of the graph to the given output stream
     */
    virtual void dump(std::ostream& out) const;

    /**
     * Add the given vertex to the graph
     * @return true if the vertex has been inserted, false otherwise
     */
    virtual bool add_vertex(uint64_t vertex_id);

    /**
     * Remove the given vertex and all edges attached to it.
     * @return true in case of success, false otherwise
     */
    virtual bool delete_vertex(uint64_t vertex_id);

    /**
     * Add the given edge in the graph
     * @return true if the edge has been inserted or updated, false in case of error
     */
    virtual bool add_edge(graph::WeightedEdge e);

    /**
     * Remove the given edge from the graph
     * @return true if the given edge has been removed, false otherwise (e.g. this edge does not exist)
     */
    virtual bool delete_edge(graph::Edge e);

    /**
     * Load the whole graph representation from the given path
     */
    virtual void load(const std::string& path);
};

} // namespace

