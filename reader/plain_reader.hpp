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
#include "reader.hpp"


namespace graph { class WeightedEdge; } // forward decl.

namespace reader {

/**
 * Specific implementation to load the edges from a plain representation (.wel). The file consists
 * of an edge per line: src dst weight.
 */
class PlainWeightedReader : public Reader {
    std::fstream m_handle; // internal file handle

public:
    /**
     * Read the edge list from the given file
     */
    PlainWeightedReader(const std::string& path);

    /**
     * Destructor
     */
    ~PlainWeightedReader();

    // interface
    bool read(graph::WeightedEdge&) override;
};


} // namespace reader


