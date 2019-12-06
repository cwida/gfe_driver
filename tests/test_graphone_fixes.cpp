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

#if defined(HAVE_GRAPHONE)

#include <cstdlib>
#include <iostream>
#include <memory>
#include <thread>
#include <vector>

#include "common/system.hpp"
#include "configuration.hpp"
#include "graph/edge.hpp"
#include "graph/edge_stream.hpp"
#if defined(HAVE_LLAMA)
#include "library/llama/llama_class.hpp"
#include "library/llama-dv/llama-dv.hpp"
#endif
#include "library/graphone/graphone.hpp"
#include "utility/graphalytics_validate.hpp"

// Log to stdout
#undef LOG
#define LOG(message) { std::cout << "\033[0;32m" << "[          ] " << "\033[0;0m" << message << std::endl; }

using namespace common::concurrency;
using namespace gfe::graph;
using namespace gfe::library;
using namespace std;

[[maybe_unused]] static std::unique_ptr<gfe::graph::WeightedEdgeStream> generate_edge_stream(uint64_t max_vector_id = 8){
    vector<gfe::graph::WeightedEdge> edges;
    for(uint64_t i = 1; i < max_vector_id; i++){
        for(uint64_t j = i + 2; j < max_vector_id; j+=2){
            edges.push_back(gfe::graph::WeightedEdge{i, j, static_cast<double>(j * 1000 + i)});
        }
    }
    return make_unique<gfe::graph::WeightedEdgeStream>(edges);
}

// Get the path to non existing temporary file
[[maybe_unused]] static string temp_file_path(){
    char pattern[] = "/tmp/gfe_XXXXXX";
    int fd = mkstemp(pattern);
    if(fd < 0){ ERROR("Cannot obtain a temporary file"); }
    close(fd); // we're going to overwrite this file anyway
    return string(pattern);
}

#if defined(HAVE_LLAMA)
/**
 * Check whether the results of GraphOne for the PageRank yield the same results of LLAMA, as they rely on the same algorithm
 */
TEST(GraphOne, PageRank) {
    LLAMAClass llama { /* directed = */ false };
    GraphOne graphone { /* directed = */ false, /* vertex dict = */ true, /* blind writes = */ true, /* num vertices */ 64 * 1024 * 1024 };
    uint64_t num_vertices = 128;

    LOG("Init stream");
    auto stream = generate_edge_stream(num_vertices);
    stream->permute();

    LOG("Insert vertices ...")
    for(uint64_t i = 1; i < num_vertices; i++){
        llama.add_vertex(i);
        graphone.add_vertex(i);
    }

    LOG("Insert edges ...");
    for(uint64_t i =0; i < stream->num_edges(); i++){
        auto edge = stream->get(i);
        llama.add_edge(edge);
        graphone.add_edge(edge);
    }

    LOG("Build ...");
    llama.build();
    graphone.build();

//    llama.dump();
//    graphone.dump();

    auto llama_results = temp_file_path();
    LOG("LLAMA PageRank: " << llama_results);
    uint64_t num_iterations = 1;
    llama.pagerank(num_iterations, 0.85, llama_results.c_str());

    auto graphone_results = temp_file_path();
    LOG("GraphOne PageRank: " << graphone_results);
    graphone.pagerank(num_iterations, 0.85, graphone_results.c_str());

    LOG("Validate the result ...");
    gfe::utility::GraphalyticsValidate::pagerank(graphone_results, llama_results);
    LOG("Validation succeeded");
}
#endif

TEST(GraphOne, ConcurrentUpdates) {
    GraphOne graphone { /* directed = */ false, /* vertex dict = */ true, /* blind writes = */ true, /* num vertices */ 64 * 1024 * 1024 };
    uint64_t num_vertices = 1<<14;
    uint64_t num_threads = 8;

    LOG("Init stream");
    auto stream = generate_edge_stream(num_vertices);
    stream->permute();

    LOG("Insert " << num_vertices << " vertices ... ");
    for(uint64_t i = 1; i < num_vertices; i++){
        graphone.add_vertex(i);
    }

    LOG("Insert " << stream->num_edges() << " edges with " << num_threads << " threads ...");
    auto worker_insert = [&graphone, &stream](uint64_t start, uint64_t end){
        for(uint64_t i = start; i < end; i++){
            auto edge = stream->get(i);
            graphone.add_edge(edge);
        }
    };

    uint64_t partition_sz = stream->num_edges() / num_threads;
    uint64_t partition_mod = stream->num_edges() % num_threads;
    uint64_t partition_start = 0;
    vector<thread> workers;
    for(uint64_t i = 0; i < num_threads; i++){
        uint64_t partition_end = partition_start + partition_sz + (i < partition_mod);
        workers.emplace_back(worker_insert, partition_start, partition_end);
        partition_start = partition_end; // next iteration
    }

    for(auto& w: workers) w.join();
    workers.clear();

    LOG("Build ...");
    graphone.build();

    // It's enough it does not end with a seg fault...
    LOG("Done");
}


#else
#include <iostream>
TEST(GraphOne, Disabled) {
    std::cout << "Tests disabled as the build does not contain the support for GraphOne." << std::endl;
}
#endif
