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

#include <cassert>

// these definitions are set by CMake when compiling the executable graphone64.
// I don't know whether all of them are actually required
#define OVER_COMMIT
#define TBB
#define PLAIN_GRAPH
#define B64
#define DEL

// Internal counters for the iterator?
//#define GRAPHONE_COUNTERS

/**
 * Do not issue warnings due to the included GraphOne library
 */
#pragma GCC diagnostic push /* works for both GCC and Clang */
#pragma GCC diagnostic ignored "-Wunused-variable" /* Do not flag unused variables */
#pragma GCC diagnostic ignored "-Woverloaded-virtual" /* Shadow methods */

// GraphOne includes
// From the directory src/
#include "graph.h" // main graph instance
#include "sgraph.h" // directed and undirected weighted graphs
#include "str.h" // dictionary external vertex id to logical id
#include "type.h" // graphone typedefs
#include "typekv.h" // first container, reserved for the "metadata"
// From the directory gview/
#include "graph_view.h" // static view

#pragma GCC diagnostic pop

// Globals (...)
extern graph* g; // defined in base.cpp

// Convenience macro, to retrieve the property graph of an instanced graph
inline static pgraph_t<lite_edge_t>* get_graphone_graph(){
   assert(g->cf_count == 2 && "Graph not initialised");
   return (pgraph_t<lite_edge_t>*) g->get_sgraph(1);
}
