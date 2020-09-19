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

#include "gtest/gtest.h"

#include <string>

#include "common/filesystem.hpp"
#include "graph/edge_stream.hpp"
#include "library/baseline/adjacency_list.hpp"
#include "library/baseline/csr.hpp"

using namespace gfe::library;
using namespace std;

// Check that all edges are properly loaded from graphs/ldbc_graphalytics/example-directed.properties
TEST(CSR, LoadDirected){
    string graph_path = common::filesystem::directory_executable() + "/graphs/ldbc_graphalytics/example-directed.properties";

    CSR csr { /* directed */ true };
    csr.load(graph_path);
    //csr.dump();

    gfe::graph::WeightedEdgeStream stream { graph_path };
    ASSERT_EQ( csr.num_edges(), stream.num_edges() );
    for(uint64_t i = 0; i < stream.num_edges(); i++){
        auto edge = stream.get(i);
        ASSERT_TRUE( csr.has_vertex(edge.source()) );
        ASSERT_TRUE( csr.has_vertex(edge.destination()) );
        ASSERT_TRUE( csr.has_edge(edge.source(), edge.destination()) );
        ASSERT_EQ( csr.get_weight(edge.source(), edge.destination()), edge.weight() );
    }
}

// Check that all edges are properly loaded from graphs/ldbc_graphalytics/example-undirected.properties
TEST(CSR, LoadUndirected){
    string graph_path = common::filesystem::directory_executable() + "/graphs/ldbc_graphalytics/example-undirected.properties";

    CSR csr { /* directed */ false };
    csr.load(graph_path);
    //csr.dump();

    gfe::graph::WeightedEdgeStream stream { graph_path };
    ASSERT_EQ( csr.num_edges(), stream.num_edges() );
    for(uint64_t i = 0; i < stream.num_edges(); i++){
        auto edge = stream.get(i);
        ASSERT_TRUE( csr.has_vertex(edge.source()) );
        ASSERT_TRUE( csr.has_vertex(edge.destination()) );
        ASSERT_TRUE( csr.has_edge(edge.source(), edge.destination()) );
        ASSERT_EQ( csr.get_weight(edge.source(), edge.destination()), edge.weight() );
        ASSERT_TRUE( csr.has_edge(edge.destination(), edge.source()) );
        ASSERT_EQ( csr.get_weight(edge.destination(), edge.source()), edge.weight() );
    }
}

