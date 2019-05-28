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

#include <string>

namespace reader {

// The format of the given graph
enum class Format {
    UNKNOWN, // No idea
    PLAIN, // Extension .el, a list with one edge per line: src dst
    PLAIN_WEIGHTED, // Extension .wel, a list with one edge per line: src dst weight
    METIS, // Extension .metis or .graph, METIS v5.1 format
};


// Retrieve the representation of the given stored graph, according to the extension of the file
Format get_graph_format(const char* path);
Format get_graph_format(const std::string& path);

} // namespace reader


