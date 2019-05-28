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

#include <memory>

namespace graph {
    class WeightedEdge; // forward decl.
}

namespace reader {

/**
 * Forward iterator to read the list of weighted edges from a stored file. Usage:
 *
 * auto reader = Reader::open( path_to_read );
 * WeightedEdge edge;
 * while ( reader.read(&edge) ) {
 *      // ... process edge ...
 * }
 *
 */
class Reader {
    // disable copy ctors
    Reader(const Reader&) = delete;
    Reader& operator=(const Reader&) = delete;

public:
    // Factory method, to obtain a specific reader depending on the extention of the given file name
    static std::unique_ptr<Reader> open(const std::string& path);

    // Default constructor
    Reader();

    // Virtual destructor
    virtual ~Reader();

    // Retrieve the next edge of the file. Returns true if an edge has been read, false if we reached the end of the file.
    bool read(graph::WeightedEdge& edge) = 0;
};

} // namespace graph
