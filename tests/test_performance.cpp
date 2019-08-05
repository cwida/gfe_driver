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

#include <cstdio>
#include <iostream>
#include <string>

#include "common/timer.hpp"
#include "experiment/insert_only.hpp"
#include "graph/edge_stream.hpp"
#include "library/baseline/adjacency_list.hpp"

using namespace common;
using namespace experiment;
using namespace std;
using namespace library;

constexpr static int num_threads = 1;
constexpr static bool is_directed = false; // whether the graph is directed
static string path_graph = "/home/dean/workspace/graphalytics/datasets/dota-league/dota-league.properties";

TEST(Performance, InsertOnly) {
    auto impl = make_shared<AdjacencyList>(is_directed);
    Timer timer;

    cout << "[Performance::InsertOnly] Loading the graph from `" << path_graph << "' ... \n";
    timer.start();
    auto stream = make_shared<graph::WeightedEdgeStream>(path_graph);
    timer.stop();
    cout << "Graph loaded in " << timer << "\n";
    stream->permute(1910);


    cout << "[Performance::InsertOnly] Executing the insertions ...\n";
    timer.start();
    InsertOnly experiment(impl, move(stream), num_threads);
    experiment.execute();
    timer.stop();
    cout << "Execution completed in " << timer << "\n";
}
