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
#include "library/interface.hpp"

namespace gfe::library {

/**
 * Wrapper to evaluate the GraphOne library
 */
class TeseoDriver : public virtual UpdateInterface {
    TeseoDriver(const TeseoDriver&) = delete;
    TeseoDriver& operator=(const TeseoDriver&) = delete;

    void* m_pImpl; // pointer to the teseo library
    const bool m_is_directed; // whether the underlying graph is directed or undirected
    std::chrono::seconds m_timeout { 0 }; // the budget to complete each of the algorithms in the Graphalytics suite

public:
    /**
     * Create an instance of Teseo
     * @param is_directed: whether the underlying graph should be directed or undirected
     */
    TeseoDriver(bool is_directed);

    /**
     * Destructor
     */
    virtual ~TeseoDriver();

    /**
     * Dump the content of the graph to given stream
     */
    void dump_ostream(std::ostream& out) const;

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
     * Returns the weight of the given edge is the edge is present, or NaN otherwise
     */
    virtual double get_weight(uint64_t source, uint64_t destination) const;

    /**
     * Check whether the graph is directed
     */
    virtual bool is_directed() const;

    /**
     * Impose a timeout on each graph computation. A computation that does not terminate by the given seconds will raise a TimeoutError.
     */
    virtual void set_timeout(uint64_t seconds);

    /**
     * Add the given vertex to the graph
     * @return true if the vertex has been inserted, false otherwise (that is, the vertex already exists)
     */
    virtual bool add_vertex(uint64_t vertex_id);

    /**
     * Remove the mapping for a given vertex. The actual internal vertex is not removed from the adjacency list.
     * @param vertex_id the vertex to remove
     * @return true if a mapping for that vertex existed, false otherwise
     */
    virtual bool remove_vertex(uint64_t vertex_id);

    /**
     * Add the given edge in the graph. The implementation does not check whether this edge already exists,
     * adding a new edge always.
     * @return always true when both the source & the destination vertices already exist, false otherwise
     */
    virtual bool add_edge(gfe::graph::WeightedEdge e);

    /**
     * Remove the given edge from the graph. There is no way to check whether the operation actually succeeded
     * in this implementation of GraphOne. Attempting to remove an edge that does not exist may result in a crash.
     * @return always true when both the source & the destination vertices already exist, false otherwise
     */
    virtual bool remove_edge(gfe::graph::Edge e);

    /**
     * Callback, invoked when a thread is created
     */
    virtual void on_thread_init(int thread_id);

    /**
     * Callback, invoked when a thread is going to be removed
     */
    virtual void on_thread_destroy(int thread_id);


};

}
