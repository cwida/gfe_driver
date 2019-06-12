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

namespace library {

// Generic exception thrown by the Stinger wrapper
DEFINE_EXCEPTION(StingerError);
    
class Stinger : public UpdateInterface, LoaderInterface {
    void* m_stinger_graph {nullptr}; // opaque object, container of the handle to the stinger graph

public:

    /**
     * Initialise the graph instance
     */
    Stinger();

    /**
     * Destructor
     */
    ~Stinger();

    /**
     * Get the number of edges contained in the graph
     */
    virtual uint64_t num_edges() const;

    /**
     * Get the number of nodes stored in the graph
     */
    virtual uint64_t num_vertices() const;

    /**
     * Returns true if the given vertex is present, false otherwise
     */
    virtual bool has_vertex(uint64_t vertex_id) const;

    /**
     * Returns true if the given edge is present, false otherwise
     */
    virtual bool has_edge(uint64_t source, uint64_t destination) const;

    /**
     * Dump the content of the graph to stdout
     */
    virtual void dump() const;

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

}

