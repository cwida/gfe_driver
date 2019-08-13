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
#include "common/error.hpp"
#include "common/filesystem.hpp"
#include "configuration.hpp"
#include "experiment/aging.hpp"
#include "graph/edge_stream.hpp"
#include "library/baseline/adjacency_list.hpp"
#include "library/stinger/stinger.hpp"
#include "library/interface.hpp"
#include "reader/graphalytics_reader.hpp"
#include "utility/graphalytics_validate.hpp"

using namespace std;
using namespace library;
using namespace utility;

const std::string path_example_directed = common::filesystem::directory_executable() + "/graphs/ldbc_graphalytics/example-directed"; // without the ext .properties at the end
const std::string path_example_undirected = common::filesystem::directory_executable() + "/graphs/ldbc_graphalytics/example-undirected"; // without the ext .properties at the end

// The algorithms to validate for a given implementation
#define GA_BFS          0x01
#define GA_PAGERANK     0x02
#define GA_WCC          0x04 /* Weakly connected components */
#define GA_CDLP         0x08 /* Community detection using Label Propagation */
#define GA_LCC          0x10 /* Local clustering coefficient */
#define GA_SSSP         0x20 /* Single source shortest paths */

// Log to stdout
#undef LOG
#define LOG(message) { std::cout << "\033[0;32m" << "[          ] " << "\033[0;0m" << message << std::endl; }

static void load_graph(library::UpdateInterface* interface, const std::string& path_graphalytics_graph){
    interface->on_main_init(1);
    interface->on_thread_init(0);

    reader::GraphalyticsReader reader { path_graphalytics_graph + ".properties" };
    uint64_t vertex_id = 0;
    while(reader.read_vertex(vertex_id)){ interface->add_vertex(vertex_id); }

    graph::WeightedEdge edge;
    while(reader.read_edge(edge)){ interface->add_edge(edge); }

    interface->on_thread_destroy(0);
    interface->on_main_destroy();
}

// Get the path to non existing temporary file
static string temp_file_path(){
    char buffer[L_tmpnam];
    tmpnam(buffer);
    return string{buffer};
}

static void validate(library::GraphalyticsInterface* interface, const std::string& path_graphalytics_graph, int algorithms = /* all */ GA_BFS | GA_PAGERANK | GA_WCC | GA_CDLP | GA_LCC | GA_SSSP){
    interface->on_main_init(1);
    interface->on_thread_init(0);

    reader::GraphalyticsReader reader { path_graphalytics_graph + ".properties" }; // to parse the properties in the file

    if(algorithms & GA_BFS){
        string path_reference = path_graphalytics_graph + "-BFS";
        if(!common::filesystem::file_exists(path_reference)) ERROR("The reference output file `" << path_reference << "' does not exist!");
        string path_result = temp_file_path();
        LOG("BFS, result: " << path_result << ", reference output: " << path_reference);
        string str_source_vertex = reader.get_property("bfs.source-vertex");
        interface->bfs(stoull(str_source_vertex, nullptr, 10), path_result.c_str());
        GraphalyticsValidate::bfs(path_result, path_reference);
        LOG("BFS, validation succeeded");
    }

    if(algorithms & GA_PAGERANK){
        string path_reference = path_graphalytics_graph + "-PR";
        if(!common::filesystem::file_exists(path_reference)) ERROR("The reference output file `" << path_reference << "' does not exist!");
        string path_result = temp_file_path();
        LOG("PAGERANK, result: " << path_result << ", reference output: " << path_reference);
        uint64_t num_iterations = stoull ( reader.get_property("pr.num-iterations") );
        double damping_factor = stod( reader.get_property("pr.damping-factor") );
        LOG("PAGERANK, parameters: num_iterations: " << num_iterations << ", damping_factor: " << damping_factor);
        interface->pagerank(num_iterations, damping_factor, path_result.c_str());
        GraphalyticsValidate::pagerank(path_result, path_reference);
        LOG("PAGERANK, validation succeeded");
    }

    if(algorithms & GA_WCC){ /* Weakly connected components */
        string path_reference = path_graphalytics_graph + "-WCC";
        if(!common::filesystem::file_exists(path_reference)) ERROR("The reference output file `" << path_reference << "' does not exist!");
        string path_result = temp_file_path();
        LOG("WCC, result: " << path_result << ", reference output: " << path_reference);
        interface->wcc(path_result.c_str());
        GraphalyticsValidate::wcc(path_result, path_reference);
        LOG("WCC, validation succeeded");
    }

    if(algorithms & GA_LCC){ /* Local clustering coefficient */
        string path_reference = path_graphalytics_graph + "-LCC";
        if(!common::filesystem::file_exists(path_reference)) ERROR("The reference output file `" << path_reference << "' does not exist!");
        string path_result = temp_file_path();
        LOG("LCC, result: " << path_result << ", reference output: " << path_reference);
        interface->lcc(path_result.c_str());
        GraphalyticsValidate::lcc(path_result, path_reference);
        LOG("LCC, validation succeeded");
    }

    if(algorithms & GA_CDLP){ /* Community detection via label propagation */
        string path_reference = path_graphalytics_graph + "-CDLP";
        if(!common::filesystem::file_exists(path_reference)) ERROR("The reference output file `" << path_reference << "' does not exist!");
        string path_result = temp_file_path();
        LOG("CDLP, result: " << path_result << ", reference output: " << path_reference);
        uint64_t max_iterations = stoull ( reader.get_property("cdlp.max-iterations") );
        LOG("CDLP, parameters: max_iterations: " << max_iterations);
        interface->cdlp(max_iterations, path_result.c_str());
        GraphalyticsValidate::cdlp(path_result, path_reference);
        LOG("CDLP, validation succeeded");
    }

    if(algorithms & GA_SSSP){ /* Dijkstra */
        string path_reference = path_graphalytics_graph + "-SSSP";
        if(!common::filesystem::file_exists(path_reference)) ERROR("The reference output file `" << path_reference << "' does not exist!");
        string path_result = temp_file_path();
        LOG("SSSP, result: " << path_result << ", reference output: " << path_reference);
        uint64_t source_vertex = stoull ( reader.get_property("sssp.source-vertex") );
        LOG("SSSP, parameters: source vertex: " << source_vertex);
        interface->sssp(source_vertex, path_result.c_str());
        GraphalyticsValidate::sssp(path_result, path_reference);
        LOG("SSSP, validation succeeded");
    }

    interface->on_thread_destroy(0);
    interface->on_main_destroy();
}

TEST(AdjacencyList, GraphalyticsDirected){
    auto adjlist = make_unique<AdjacencyList>(/* directed */ true);
    load_graph(adjlist.get(), path_example_directed);
    validate(adjlist.get(), path_example_directed);
}

TEST(AdjacencyList, GraphalyticsUndirected){
    auto adjlist = make_unique<AdjacencyList>(/* directed */ false);
    load_graph(adjlist.get(), path_example_undirected);
    validate(adjlist.get(), path_example_undirected);
}

TEST(AdjacencyList, Aging){
    auto adjlist = make_shared<AdjacencyList>(/* directed */ false);
    auto edge_stream = make_shared<graph::WeightedEdgeStream>(path_example_undirected + ".properties" );
    experiment::Aging aging(adjlist, move(edge_stream), 1024 * 32, 8);
    aging.execute();
    validate(adjlist.get(), path_example_undirected);
}

#if defined(HAVE_STINGER)
TEST(Stinger, GraphalyticsDirected){
    auto graph = make_unique<Stinger>(/* directed */ true);
    load_graph(graph.get(), path_example_directed);
    validate(graph.get(), path_example_directed);
}

TEST(Stinger, GraphalyticsUndirected){
    auto graph = make_unique<Stinger>(/* directed */ false);
    load_graph(graph.get(), path_example_undirected);
    validate(graph.get(), path_example_undirected);
}
#endif
