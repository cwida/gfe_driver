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
#include "experiment/aging2_experiment.hpp"
#include "graph/edge_stream.hpp"
#include "library/baseline/adjacency_list.hpp"

using namespace experiment;
using namespace graph;
using namespace library;
using namespace std;

static
void validate_aging2(bool is_directed, const string& path_graph, const string& path_log, uint64_t exp_granularity = 1024){
    auto stream = make_shared<WeightedEdgeStream>(path_graph);
    auto adjlist = make_shared<AdjacencyList>(is_directed);

    Aging2Experiment exp_aging;
    exp_aging.set_library(adjlist);
    exp_aging.set_log(path_log);
    exp_aging.set_parallelism_degree(8);
    exp_aging.set_worker_granularity(exp_granularity);
    exp_aging.execute();

    adjlist->dump();

    stream = make_shared<WeightedEdgeStream>(path_graph);
    unordered_set<uint64_t> vertices;
    for(uint64_t i = 0, sz = stream->num_edges(); i < sz; i++){
        auto edge = stream->get(i);
        ASSERT_TRUE(adjlist->has_vertex(edge.source()));
        ASSERT_TRUE(adjlist->has_vertex(edge.destination()));
//        cout << "check edge: " << edge.source() << ", destination: " << edge.destination() << ", weight: " << edge.weight() << endl;
        ASSERT_TRUE(adjlist->has_edge(edge.source(), edge.destination()));
        double weight = adjlist->get_weight(edge.source(), edge.destination());
        ASSERT_DOUBLE_EQ(weight, edge.weight());

        vertices.insert(edge.source());
        vertices.insert(edge.destination());
    }

    ASSERT_EQ(stream->num_edges(), adjlist->num_edges());
    ASSERT_EQ(vertices.size(), adjlist->num_vertices());
}


TEST(Aging2, Undirected){
    const string path_graph = common::filesystem::directory_executable() + "/graphs/ldbc_graphalytics/example-undirected.properties";
    const string path_log = common::filesystem::directory_executable() + "/graphs/ldbc_graphalytics/example-undirected.graphlog";
    validate_aging2(/* is directed ? */ false, path_graph, path_log, 4);
}
