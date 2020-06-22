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

#if defined(HAVE_LLAMA)
#include "library/llama/llama_internal.hpp"
#include "library/llama/llama_ref.hpp"
#endif

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

static const uint64_t num_iterations_default = 1ull << 10;
[[maybe_unused]] static uint64_t get_num_iterations(){
    static const char* gfe_num_iterations =  getenv("GFE_NUM_ITERATIONS");
    if(gfe_num_iterations == nullptr){
        cout << "Warning: env. var. GFE_NUM_ITERATIONS not set. Using the default: " << num_iterations_default << endl;
        return num_iterations_default;
    } else {
        return strtoull(gfe_num_iterations, nullptr, 10);
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

#if defined(HAVE_LLAMA)
static void _test_perf_run_llama(){
    auto impl = make_shared<LLAMAClass>(is_directed);
    Timer timer;
    string path_graph = get_path_graph();

    cout << "[Performance::LLAMA_Iterator_Overhead] Loading the graph from `" << path_graph << "' ... \n";
    timer.start();
    auto stream = make_shared<gfe::graph::WeightedEdgeStream>(path_graph);
    timer.stop();
    cout << "Graph loaded in " << timer << "\n";
    stream->permute(1910);


    cout << "[Performance::LLAMA_Iterator_Overhead] Executing the insertions ...\n";
    timer.start();
    InsertOnly experiment(impl, move(stream), num_threads);
    experiment.execute();
    timer.stop();
    cout << "Execution completed in " << timer << "\n";


    uint64_t num_iterations = get_num_iterations();
    uint64_t num_vertices = impl->num_vertices();
    cout << "[Performance::LLAMA_Iterator_Overhead] Initialising " << num_iterations << " the output iterator...\n";
    auto instance = dynamic_cast<LLAMAClass*>(impl.get());
    shared_lock<LLAMAClass::shared_mutex_t> slock(instance->m_lock_checkpoint);
    auto graph = instance->get_snapshot();
    timer.start();
    for(uint64_t i = 0; i < num_iterations; i++){
        ll_edge_iterator iterator;
        graph.out_iter_begin(iterator, i % num_vertices);
    }
    timer.stop();
    double time_per_item_usecs = static_cast<double>(timer.microseconds()) / num_iterations;
    cout << "Execution completed in " << timer << ". Init time per item: " << time_per_item_usecs << " microseconds\n";

    cout << "[Performance::LLAMA_Iterator_Overhead] Initialising " << num_iterations << " the input iterator...\n";
    timer.start();
    for(uint64_t i = 0; i < num_iterations; i++){
        ll_edge_iterator iterator;
        graph.out_iter_begin(iterator, i % num_vertices);
    }
    timer.stop();
    time_per_item_usecs = static_cast<double>(timer.microseconds()) / num_iterations;
    cout << "Execution completed in " << timer << ". Init time per item: " << time_per_item_usecs << " microseconds\n";
}

TEST(Performance, LLAMA_Iterator_Overhead) {
    _test_perf_run_llama();
}
#endif
