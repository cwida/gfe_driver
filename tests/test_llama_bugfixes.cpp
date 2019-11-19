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

#include "common/system.hpp"
#include "configuration.hpp"
#include "graph/edge.hpp"
#include "library/interface.hpp"
#include "library/llama/llama_class.hpp"
#include "library/llama-dv/llama-dv.hpp"

// Log to stdout
#undef LOG
#define LOG(message) { std::cout << "\033[0;32m" << "[          ] " << "\033[0;0m" << message << std::endl; }

using namespace common::concurrency;
using namespace graph;
using namespace library;
using namespace std;

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

TEST(LLAMA, OutDegreeCounterOverflow) { // This is going to take a while, ZzZ...
    LLAMAClass llama { /* directed = */ false };
    llama.add_vertex(1);

    uint64_t num_vertices = 1ull << 17; // more than 65k
    uint64_t num_threads = 8;

    LOG("Insertions ...");
    vector<thread> threads;
    for(uint64_t thread_id = 0; thread_id < num_threads; thread_id++){
        threads.emplace_back([num_threads, num_vertices, &llama](uint64_t thread_id){
            set_thread_name("Worker #" + to_string(thread_id));

            for(uint64_t j = 2 + thread_id; j < num_vertices; j += num_threads){
                llama.add_vertex(j);
                WeightedEdge edge{ 1, j, 10.0 * j};

                ASSERT_TRUE ( llama.add_edge( edge ) );
            }

        }, thread_id);
    }
    for(auto& t: threads) t.join();
    threads.clear();

    LOG("Updates (% 101)..."); // this block is not necessary for the bug, just to make the test more complex
    for(uint64_t thread_id = 0; thread_id < num_threads; thread_id++){
        threads.emplace_back([num_threads, num_vertices, &llama](uint64_t thread_id){
            set_thread_name("Worker #" + to_string(thread_id));

            for(uint64_t j = 2 + thread_id; j < num_vertices; j += num_threads){
                if(j % 101 != 0) continue;
                Edge edge_old {1, j};
                ASSERT_TRUE( llama.remove_edge(edge_old) );

                WeightedEdge edge_new{ 1, j, 100.0 * j};
                ASSERT_TRUE ( llama.add_edge( edge_new ) );
            }

        }, thread_id);
    }
    for(auto& t: threads) t.join();
    threads.clear();

    LOG("Build ...");
    llama.build(); // <- CRASH

    LOG("Updates ...");
    for(uint64_t thread_id = 0; thread_id < num_threads; thread_id++){
        threads.emplace_back([num_threads, num_vertices, &llama](uint64_t thread_id){
            set_thread_name("Worker #" + to_string(thread_id));

            for(uint64_t j = 2 + thread_id; j < num_vertices; j += num_threads){
                Edge edge_old {1, j};
                ASSERT_TRUE( llama.remove_edge(edge_old) );

                WeightedEdge edge_new{ 1, j, 1000.0 * j};
                ASSERT_TRUE ( llama.add_edge( edge_new ) );
            }

        }, thread_id);
    }

    for(auto& t: threads) t.join();
    threads.clear();
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
