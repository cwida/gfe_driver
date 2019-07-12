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
 * Loader for the DIMACS challenge #9 file format, year 2005/06, spec: http://users.diag.uniroma1.it/challenge9/format.shtml#graph
 *
 * The file always represents a weighted graph. The format is ASCII and it is similar to plain edge list. The first character of
 * each line determines the type of information that follows on the same line:
 * c - the rest of the line is a comment and should be ignored
 * p - problem specification, ignored for this loader
 * a - an edge in the format src dst weight
 *
 * Example, a graph with three edges 1->2 (weight=10), 1->5 (weight=100), 2->3 (weight=20):
 * c This is a header comment
 * c and should be ignored
 * p sp n m
 * a 1 2 10
 * a 1 5 100
 * a 2 3 20
 */
class Dimacs9Reader : public Reader {
    std::fstream m_handle; // internal file handle
    char m_buffer[1025]; // read buffer

    /**
     * Check whether the given line should be parsed, that is, it represents an edge.
     * @param line the content to check. As side effect, the pointer is moved to the next character
     */
    static bool parse_line(char*& line);

    /**
     * Check whether the current marker points to a number
     */
    static bool is_number(const char* marker);
public:
    /**
     * Read the edge list from the given file
     * @param path the source of the file
     */
    Dimacs9Reader(const std::string& path);

    /**
     * Destructor
     */
    ~Dimacs9Reader();

    /**
     * Interface, report one edge at the time
     */
    bool read(graph::WeightedEdge&) override;

    /**
     * Interface, whether the input graph is directed
     */
    bool is_directed() const override;
};


} // namespace
