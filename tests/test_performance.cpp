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
#include <cstdlib> // getenv
#include <iostream>
#include <string>

#include "common/system.hpp"
#include "common/timer.hpp"
#include "experiment/insert_only.hpp"
#include "graph/edge_stream.hpp"
#include "library/baseline/adjacency_list.hpp"

using namespace common;
using namespace gfe::experiment;
using namespace gfe::library;
using namespace std;

extern char** environ;

constexpr static int num_threads = 1;
constexpr static bool is_directed = false; // whether the graph is directed

static const string path_graph_default = common::filesystem::directory_executable() + "/graphs/ldbc_graphalytics/example-undirected.properties";
static string get_path_graph(){
    static const char* path_graph = getenv("GFE_PATH_GRAPH");
    if(path_graph == nullptr){
        cout << "Warning: env. var. GFE_PATH_GRAPH not set. Using the default: " << path_graph_default << endl;
        return path_graph_default;
    } else {
        return string(path_graph);
    }
}


TEST(Performance, InsertOnly) {
    auto impl = make_shared<AdjacencyList>(is_directed);
    Timer timer;
    string path_graph = get_path_graph();

    cout << "[Performance::InsertOnly] Loading the graph from `" << path_graph << "' ... \n";
    timer.start();
    auto stream = make_shared<gfe::graph::WeightedEdgeStream>(path_graph);
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


TEST(Performance, LCC) {
    auto impl = make_shared<AdjacencyList>(is_directed);
    Timer timer;
    string path_graph = get_path_graph();

    cout << "[Performance::LCC] Loading the graph from `" << path_graph << "' ... \n";
    timer.start();
    auto stream = make_shared<gfe::graph::WeightedEdgeStream>(path_graph);
    timer.stop();
    cout << "Graph loaded in " << timer << "\n";
    stream->permute(1910);


    cout << "[Performance::LCC] Executing the insertions ...\n";
    timer.start();
    InsertOnly experiment(impl, move(stream), num_threads);
    experiment.execute();
    timer.stop();
    cout << "Execution completed in " << timer << "\n";

    cout << "[Performance::LCC] Executing the LCC algorithm ...\n";
    timer.start();
    impl->lcc();
    timer.stop();
    cout << "Execution completed in " << timer << "\n";

}
