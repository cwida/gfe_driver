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
#include <cstdlib> // mkstemp
#include <iostream>
#include <random>
#include <string>
#include <unordered_set>

#include "common/error.hpp"
#include "common/filesystem.hpp"
#include "common/permutation.hpp"
#include "configuration.hpp"
#include "graph/edge_stream.hpp"
#include "library/baseline/adjacency_list.hpp"
#include "library/baseline/csr.hpp"
#if defined(HAVE_LLAMA)
#include "library/llama/llama_class.hpp"
#include "library/llama/llama_ref.hpp"
#endif
#if defined(HAVE_STINGER)
#include "library/stinger/stinger.hpp"
#endif
#if defined(HAVE_GRAPHONE)
#include "library/graphone/graphone.hpp"
#endif
#if defined(HAVE_LIVEGRAPH)
#include "library/livegraph/livegraph_driver.hpp"
#endif
#if defined(HAVE_TESEO)
#include "library/teseo/teseo_driver.hpp"
#endif
#include "library/interface.hpp"
#include "reader/graphalytics_reader.hpp"
#include "utility/graphalytics_validate.hpp"

using namespace gfe::library;
using namespace gfe::utility;
using namespace std;

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

static void load_graph(gfe::library::UpdateInterface* interface, const std::string& path_graphalytics_graph){
    interface->on_main_init(1);
    interface->on_thread_init(0);

    gfe::reader::GraphalyticsReader reader { path_graphalytics_graph + ".properties" };
    uint64_t vertex_id = 0;
    while(reader.read_vertex(vertex_id)){ interface->add_vertex(vertex_id); }

    gfe::graph::WeightedEdge edge;
    while(reader.read_edge(edge)){ interface->add_edge(edge); }

    interface->build();
    interface->on_thread_destroy(0);
    interface->on_main_destroy();
}

// Get the path to non existing temporary file
static string temp_file_path(){
    char pattern[] = "/tmp/gfe_XXXXXX";
    int fd = mkstemp(pattern);
    if(fd < 0){ ERROR("Cannot obtain a temporary file"); }
    close(fd); // we're going to overwrite this file anyway
    return string(pattern);
}

static void validate(gfe::library::GraphalyticsInterface* interface, const std::string& path_graphalytics_graph, int algorithms = /* all */ GA_BFS | GA_PAGERANK | GA_WCC | GA_CDLP | GA_LCC | GA_SSSP){
    interface->on_main_init(1);
    interface->on_thread_init(0);

    gfe::reader::GraphalyticsReader reader { path_graphalytics_graph + ".properties" }; // to parse the properties in the file

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

TEST(CSR, GraphalyticsDirected){
    auto csr = make_unique<CSR>(/* directed */ true);
    csr->load(path_example_directed + ".properties");
    validate(csr.get(), path_example_directed);
}

TEST(CSR, GraphalyticsUndirected){
    auto csr = make_unique<CSR>(/* directed */ false);
    csr->load(path_example_undirected + ".properties");
    validate(csr.get(), path_example_undirected);
}

#if defined(HAVE_LLAMA)
TEST(LLAMA, GraphalyticsDirected){
    auto graph = make_unique<LLAMAClass>(/* directed */ true);
    load_graph(graph.get(), path_example_directed);
    validate(graph.get(), path_example_directed);
}

TEST(LLAMA, GraphalyticsUndirected){
    auto graph = make_unique<LLAMAClass>(/* directed */ false);
    load_graph(graph.get(), path_example_undirected);
    validate(graph.get(), path_example_undirected);
}

TEST(LLAMARef, GraphalyticsDirected){
    auto graph = make_unique<LLAMARef>(/* directed */ true);
    load_graph(graph.get(), path_example_directed);
    validate(graph.get(), path_example_directed);
}

TEST(LLAMARef, GraphalyticsUndirected){
    auto graph = make_unique<LLAMARef>(/* directed */ false);
    load_graph(graph.get(), path_example_undirected);
    validate(graph.get(), path_example_undirected);
}

#if defined(HAVE_STINGER)
TEST(LLAMARef, GraphalyticsCompareWCCWithStinger){
    auto stinger_ref = make_unique<StingerRef>(/* directed */ false);
    auto llama_ref = make_unique<LLAMARef>(/* directed */ false);

    // generate a random vertex stream
    const uint64_t max_vertex_id = 1600000;
    vector<gfe::graph::WeightedEdge> edges;
    for(uint64_t i = 10; i < max_vertex_id; i+=10){
        for(uint64_t j = i +10; j < std::min(max_vertex_id, i + 50); j+=10){
            edges.push_back(gfe::graph::WeightedEdge{i, j, static_cast<double>(j * 1000 + i)});
        }
    }
    gfe::graph::WeightedEdgeStream stream(edges);
    stream.permute();

    unordered_set<uint64_t> existing_vertices;
    auto insert_vertex = [&existing_vertices, &llama_ref, &stinger_ref](uint64_t vertex){
        if(existing_vertices.count(vertex) == 0){
            llama_ref->add_vertex(vertex);
            stinger_ref->add_vertex(vertex);
            existing_vertices.insert(vertex);
        }
    };

    LOG("Insert " << stream.num_edges() << " edges ...");
    for(uint64_t i = 0; i < stream.num_edges(); i++){
        insert_vertex(stream[i].m_source);
        insert_vertex(stream[i].m_destination);

        llama_ref->add_edge(stream[i]);
        stinger_ref->add_edge(stream[i]);

        if(i == 1000) llama_ref->build(); // create multiple snapshots
    }

    llama_ref->build();

    ASSERT_EQ(llama_ref->num_vertices(), stinger_ref->num_vertices());
    ASSERT_EQ(llama_ref->num_edges(), stinger_ref->num_edges());

    string path_result_stinger_ref = temp_file_path();
    LOG("WCC StingerRef: " << path_result_stinger_ref);
    stinger_ref->wcc(path_result_stinger_ref.c_str());

    string path_result_llama_ref = temp_file_path();
    LOG("WCC LLAMARef: " << path_result_llama_ref);
    llama_ref->wcc(path_result_llama_ref.c_str());

    LOG("Validate the result ...");
    gfe::utility::GraphalyticsValidate::wcc(path_result_llama_ref, path_result_stinger_ref);
    LOG("Validation succeeded");
}

#endif

#endif

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

TEST(StingerRef, GraphalyticsDirected){
    auto graph = make_unique<StingerRef>(/* directed */ true);
    load_graph(graph.get(), path_example_directed);
    validate(graph.get(), path_example_directed);
}

TEST(StingerRef, GraphalyticsUndirected){
    auto graph = make_unique<StingerRef>(/* directed */ false);
    load_graph(graph.get(), path_example_undirected);
    validate(graph.get(), path_example_undirected);
}

TEST(StingerRef, GraphalyticsCompareWithBaseline){
    auto stinger_ref = make_unique<StingerRef>(/* directed */ false);
    auto baseline = make_unique<AdjacencyList>(/* directed */ false);

    // generate a random vertex stream
    const uint64_t max_vertex_id = 1600;
    vector<gfe::graph::WeightedEdge> edges;
    for(uint64_t i = 10; i < max_vertex_id; i+=10){
        for(uint64_t j = i +10; j < max_vertex_id; j+=10){
            edges.push_back(gfe::graph::WeightedEdge{i, j, static_cast<double>(j * 1000 + i)});
        }
    }
    gfe::graph::WeightedEdgeStream stream(edges);
    stream.permute();

    // insert the vertices in the libraries, one by one
    LOG("Insert " << max_vertex_id / 10 << " vertices ...");
    vector<int64_t> vertices;
    for(uint64_t i = 10; i <= max_vertex_id; i+= 10){ vertices.push_back(i); }
    unique_ptr<uint64_t[]> ptr_permutation { new uint64_t[vertices.size()] };
    common::permute(ptr_permutation.get(), vertices.size(), 3);
    for(uint64_t i = 0; i < vertices.size(); i++){
        uint64_t vertex = vertices[ ptr_permutation[i] ];
        baseline->add_vertex(vertex);
        stinger_ref->add_vertex(vertex);
    }

    LOG("Insert " << stream.num_edges() << " edges ...");
    for(uint64_t i = 0; i < stream.num_edges(); i++){
        baseline->add_edge(stream[i]);
        stinger_ref->add_edge(stream[i]);
    }

    ASSERT_EQ(baseline->num_vertices(), stinger_ref->num_vertices());
    ASSERT_EQ(baseline->num_edges(), stinger_ref->num_edges());

    uint64_t source_vertex = 10;

    string path_result_baseline = temp_file_path();
    LOG("BFS Baseline: " << path_result_baseline);
    baseline->bfs(source_vertex, path_result_baseline.c_str());

    string path_result_stinger_ref = temp_file_path();
    LOG("BFS StingerRef: " << path_result_stinger_ref);
    stinger_ref->bfs(source_vertex, path_result_stinger_ref.c_str());

    LOG("Validate the result ...");
    gfe::utility::GraphalyticsValidate::bfs(path_result_stinger_ref, path_result_baseline);
    LOG("Validation succeeded");
}

#endif

#if defined(HAVE_GRAPHONE)
TEST(GraphOne, GraphalyticsDirected){
    auto graph = make_unique<GraphOne>(/* directed */ true, /* vertex dictionary */ true, /* blind writes */ true, /* ignore build = */ false, /* ref impl ? */ false, /* max num vertices */ 32);
    load_graph(graph.get(), path_example_directed);
    validate(graph.get(), path_example_directed);
}

TEST(GraphOne, GraphalyticsUndirected){
    auto graph = make_unique<GraphOne>(/* directed */ false, /* vertex dictionary */ true, /* blind writes */ true, /* ignore build = */ false, /* ref impl ? */ false, /* max num vertices */ 32);
    load_graph(graph.get(), path_example_undirected);
    validate(graph.get(), path_example_undirected);
}

TEST(GraphOneRef, GraphalyticsUndirected){
    auto graph = make_unique<GraphOne>(/* directed */ false, /* vertex dictionary */ true, /* blind writes */ true, /* ignore build = */ false, /* ref impl ? */ true, /* max num vertices */ 32);
    load_graph(graph.get(), path_example_undirected);
    validate(graph.get(), path_example_undirected);
}
#endif

#if defined(HAVE_LIVEGRAPH)
TEST(LiveGraph, GraphalyticsUndirected){
    auto graph = make_unique<LiveGraphDriver>(/* directed */ false);
    load_graph(graph.get(), path_example_undirected);
    validate(graph.get(), path_example_undirected);
}
TEST(LiveGraph, GraphalyticsDirected){
    auto graph = make_unique<LiveGraphDriver>(/* directed */ true);
    load_graph(graph.get(), path_example_directed);

    /**
     * On directed graphs, for LiveGraph we only have an implementation of WCC & SSSP. All the other kernel implementations
     * rely on the knowledge of incoming edges, which we are not keeping atm.
     */
    validate(graph.get(), path_example_directed, GA_WCC | GA_SSSP);
}
#endif

#if defined(HAVE_TESEO)
TEST(Teseo, GraphalyticsUndirected){
    auto graph = make_unique<TeseoDriver>(/* directed */ false);
    load_graph(graph.get(), path_example_undirected);
    validate(graph.get(), path_example_undirected);
}

TEST(TeseoLCC, GraphalyticsUndirected){
    auto graph = make_unique<TeseoDriverLCC>(/* directed */ false);
    load_graph(graph.get(), path_example_undirected);
    validate(graph.get(), path_example_undirected, GA_LCC);
}

#endif
