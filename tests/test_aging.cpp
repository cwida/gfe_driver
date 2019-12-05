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

#include <cstdlib> // getenv
#include <memory>
#include <unordered_set>

#include "common/filesystem.hpp"
#include "experiment/aging.hpp"
#include "graph/edge_stream.hpp"
#include "library/baseline/adjacency_list.hpp"

using namespace gfe::experiment;
using namespace gfe::graph;
using namespace gfe::library;
using namespace std;

static
void validate_aging(bool is_directed, const string& path_graph, uint64_t exp_granularity = 1024){
    auto stream = make_shared<WeightedEdgeStream>(path_graph);
    auto adjlist = make_shared<AdjacencyList>(is_directed);

    Aging exp_aging{adjlist, move(stream), /* aging coeff. */ 8.0, /* num threads */  8};
    exp_aging.set_operation_granularity(exp_granularity);
    exp_aging.execute();

    stream = make_shared<WeightedEdgeStream>(path_graph);
    unordered_set<uint64_t> vertices;
    for(uint64_t i = 0, sz = stream->num_edges(); i < sz; i++){
        auto edge = stream->get(i);
        ASSERT_TRUE(adjlist->has_vertex(edge.source()));
        ASSERT_TRUE(adjlist->has_vertex(edge.destination()));
        double weight = adjlist->get_weight(edge.source(), edge.destination());
        ASSERT_DOUBLE_EQ(weight, edge.weight());

        vertices.insert(edge.source());
        vertices.insert(edge.destination());
    }

    ASSERT_EQ(stream->num_edges(), adjlist->num_edges());
    ASSERT_EQ(vertices.size(), adjlist->num_vertices());
}

TEST(Aging, Directed){
    validate_aging(true, common::filesystem::directory_executable() + "/graphs/rome99.no_duplicates.gr", 64);
}

TEST(Aging, Undirected){
    validate_aging(false, common::filesystem::directory_executable() + "/graphs/ldbc_graphalytics/example-undirected.properties", 4);
}
