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


#include "library/interface.hpp"

namespace gfe::library {

/**
 * This is a dummy instance of the interface. All operations are *nop*. It is used simply
 * to compute the driver and network overhead in running the experiments.
 */
class Dummy : public virtual UpdateInterface {
private:
    const bool m_is_directed; // whether the underlying graph is directed
public:

    /**
     * Initialise the instance
     */
    Dummy(bool is_directed);

    /**
     * Destructor
     */
    ~Dummy();

    virtual uint64_t num_edges() const; // returns 0
    virtual uint64_t num_vertices() const; // returns 0
    virtual bool is_directed() const; // returns m_directed
    virtual bool has_vertex(uint64_t vertex_id) const; // always returns true
    virtual double get_weight(uint64_t source, uint64_t destination) const; // always returns 0
    virtual void dump_ostream(std::ostream& out) const; // prints `DuMMy'
    virtual bool add_vertex(uint64_t vertex_id); // returns true
    virtual bool remove_vertex(uint64_t vertex_id); // returns true
    virtual bool add_edge(graph::WeightedEdge e); // returns true
    virtual bool add_edge_v2(gfe::graph::WeightedEdge e); // returns true
    virtual bool remove_edge(graph::Edge e); // returns true
    virtual void set_timeout(uint64_t seconds); // nop
}; // class

} // namespace
