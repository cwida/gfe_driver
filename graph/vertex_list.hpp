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

namespace gfe::graph {

class CByteArray;

/**
 * Represent a list of vertices in a graph
 */
class VertexList {
private:
    VertexList(const VertexList&) = delete;
    VertexList& operator=(const VertexList&) = delete;

    CByteArray* m_vertices; // the list of vertices contained
public:
    /**
     * Pass the list of vertices
     */
    VertexList(CByteArray* vertices);

    /**
     * Destructor
     */
    ~VertexList();

    /**
     * Retrieve the vertex at the given position, in [0, num_vertices() )
     */
    uint64_t get(uint64_t index) const;
    uint64_t operator[](uint64_t index) const { return get(index); } // alias

    /**
     * Total number of vertices
     */
    uint64_t num_vertices() const noexcept;

    /**
     * Perform a random permutation of the list of vertices, using a custom seed
     */
    void permute();

    /**
     * Perform a random permutation of the list of vertices, using the given seed
     */
    void permute(uint64_t seed);

};
    
} // namespace
