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
#include <functional> // std::hash
#include <ostream>
#include <utility> // std::swap

namespace graph {

// Simple  representation of an edge as a pair <source, destination>
class Edge {
public:
    uint64_t m_source;
    uint64_t m_destination;

    Edge(): m_source{0}, m_destination{0} {}
    Edge(uint64_t source, uint64_t destination) : m_source(source), m_destination(destination) { };

    uint64_t source() const { return m_source; }
    uint64_t destination() const { return m_destination; }

    // Check whether the two edges are equal
    bool operator==(const Edge&) const noexcept;
    bool operator!=(const Edge&) const noexcept;

    // Swap source & destination ID for the current edge
    void swap_src_dst() noexcept { std::swap(m_source, m_destination); }
};

// Simple representation of an edge as a tuple <source, destination, weight>
class WeightedEdge : public Edge{
public:
    WeightedEdge();
    WeightedEdge(uint64_t source, uint64_t destination, double weight);

    double m_weight;
    double weight() const { return m_weight; }

    // Get a copy of the given (non weighted) edge
    Edge edge() const { return Edge{source(), destination()}; }

    // Check whether the two edges are equal
    bool operator==(const WeightedEdge&) const noexcept;
    bool operator!=(const WeightedEdge&) const noexcept;
};

std::ostream& operator<<(std::ostream& out, const Edge& e);
std::ostream& operator<<(std::ostream& out, const WeightedEdge& e);
} // namespace graph

namespace std {
template<> struct hash<graph::Edge>{ // hash function for graph::Edge
private:
    // Adapted from the General Purpose Hash Function Algorithms Library
    // Author: Arash Partow - 2002
    // URL: http://www.partow.net
    // URL: http://www.partow.net/programming/hashfunctions/index.html
    // MIT License
    unsigned int APHash(uint64_t value) const {
       unsigned int hash = 0xAAAAAAAA;
       unsigned int i    = 0;
       char* str = (char*) &value;

       for (i = 0; i < sizeof(value); ++str, ++i) {
          hash ^= ((i & 1) == 0) ? (  (hash <<  7) ^ (*str) * (hash >> 3)) :
                                   (~((hash << 11) + ((*str) ^ (hash >> 5))));
       }

       return hash;
    }

public:
    size_t operator()(const graph::Edge& e) const {
        return APHash( (e.source() << 32) + (e.destination() & 0xFFFFFFFF) );
    }
};

} // namespace std
