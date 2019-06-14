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
#include <ostream>

namespace graph {

// Simple  representation of an edge as a pair <source, destination>
class Edge {
public:
    uint64_t m_source;
    uint64_t m_destination;

    uint64_t source() const { return m_source; }
    uint64_t destination() const { return m_destination; }

    // Check whether the two edges are equal
    bool operator==(const Edge&) const noexcept;
    bool operator!=(const Edge&) const noexcept;
};

// Simple representation of an edge as a tuple <source, destination, weight>
class WeightedEdge : public Edge{
public:
    WeightedEdge();
    WeightedEdge(uint64_t source, uint64_t destination, uint32_t weight);

    uint32_t m_weight;
    uint32_t weight() const { return m_weight; }

    // Get a copy of the given (non weighted) edge
    Edge edge() const { return Edge{source(), destination()}; }

    // Check whether the two edges are equal
    bool operator==(const WeightedEdge&) const noexcept;
    bool operator!=(const WeightedEdge&) const noexcept;
};

std::ostream& operator<<(std::ostream& out, const Edge& e);
std::ostream& operator<<(std::ostream& out, const WeightedEdge& e);

} // namespace graph



