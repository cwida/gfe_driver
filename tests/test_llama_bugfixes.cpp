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

// Tests specific to the LLAMA implementation
#include "gtest/gtest.h"

#if defined(HAVE_LLAMA)

#include <iostream>
#include <memory>
#include <thread>
#include <unordered_set>
#include <vector>

#include "configuration.hpp"
#include "graph/edge.hpp"
#include "library/interface.hpp"
#include "library/llama/llama_class.hpp"
#include "library/llama-dv/llama-dv.hpp"

using namespace std;
using namespace graph;
using namespace library;

TEST(LLAMA, RemoveInsertEdge) {
    LLAMAClass llama { /* directed = */ false };

    ASSERT_TRUE( llama.add_vertex(10) );
    ASSERT_TRUE( llama.add_vertex(20) );
    ASSERT_TRUE( llama.add_vertex(30) );
    ASSERT_TRUE( llama.add_edge(WeightedEdge{10, 20, 1020}));
    ASSERT_TRUE( llama.has_edge(10, 20 ) );
    ASSERT_TRUE( llama.has_edge(20, 10 ) ); // because the graph is undirected

    llama.build();

    ASSERT_TRUE( llama.has_edge(10, 20 ) );
    ASSERT_TRUE( llama.has_edge(20, 10 ) ); // because the graph is undirected
    ASSERT_EQ ( (int) llama.get_weight(10, 20), 1020 );

    ASSERT_TRUE( llama.remove_edge(Edge{10, 20}));
    ASSERT_FALSE( llama.has_edge(10, 20 ) );
    ASSERT_FALSE( llama.has_edge(20, 10 ) ); // because the graph is undirected

    ASSERT_TRUE( llama.add_edge(WeightedEdge{20, 10, 2010}) );
    ASSERT_FALSE(llama.add_edge(WeightedEdge{10, 20, 10201}) );

    ASSERT_TRUE( llama.has_edge(10, 20 ) );
    ASSERT_TRUE( llama.has_edge(20, 10 ) ); // because the graph is undirected
    ASSERT_EQ ( (int) llama.get_weight(10, 20), 2010 );
    ASSERT_EQ ( (int) llama.get_weight(20, 10), 2010 );

    llama.build();
    ASSERT_TRUE( llama.has_edge(10, 20 ) );
    ASSERT_TRUE( llama.has_edge(20, 10 ) );
    ASSERT_EQ ( (int) llama.get_weight(10, 20), 2010 );
    ASSERT_EQ ( (int) llama.get_weight(20, 10), 2010 );
}

TEST(LLAMA, RewriteEdgeInWriteStore) {
    LLAMAClass llama { /* directed = */ false };

    ASSERT_TRUE( llama.add_vertex(10) );
    ASSERT_TRUE( llama.add_vertex(20) );
    ASSERT_TRUE( llama.add_vertex(30) );
    ASSERT_TRUE( llama.add_edge(WeightedEdge{10, 20, 1020}));

    llama.build();

    ASSERT_FALSE( llama.remove_edge(Edge{20, 30}) );

    ASSERT_TRUE( llama.add_edge(WeightedEdge{20, 30, 2030}) );
    ASSERT_TRUE( llama.has_edge(20, 30 ) );
//    ASSERT_EQ( (int) llama.get_weight(20, 30), 2030 ); // only visible after build?

    ASSERT_TRUE( llama.remove_edge(Edge{20, 30}) );
    ASSERT_FALSE(  llama.has_edge(20, 30 ) );
    ASSERT_FALSE(  llama.has_edge(30, 20 ) );
    ASSERT_FALSE( llama.remove_edge(Edge{20, 30}) ); // already removed

    ASSERT_TRUE( llama.add_edge(WeightedEdge{30, 20, 3020}) ); // undirected
    ASSERT_TRUE( llama.has_edge(30, 20) );
//    ASSERT_EQ( (int) llama.get_weight(30, 20), 3020 ); // only visible after build?

    llama.build();

    ASSERT_TRUE( llama.has_edge(30, 20) );
    ASSERT_EQ( (int) llama.get_weight(30, 20), 3020 );
    ASSERT_EQ( (int) llama.get_weight(20, 30), 3020 );
}

TEST(LLAMA_DV, RemoveEdge) {
    LLAMA_DV llama { /* directed = */ false };

    ASSERT_TRUE( llama.add_vertex(1) );
    ASSERT_TRUE( llama.add_vertex(2) );
    ASSERT_TRUE( llama.add_vertex(3) );
    ASSERT_TRUE( llama.add_edge(WeightedEdge{1, 2, 1020}));
    ASSERT_TRUE( llama.remove_edge(Edge{1, 2}));

}

#else
#include <iostream>
TEST(LLAMA, Disabled) {
    std::cout << "Tests disabled as the build does not contain the support for LLAMA." << std::endl;
}
#endif
