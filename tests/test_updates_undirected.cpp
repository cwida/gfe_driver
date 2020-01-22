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

/**
 * Test the update interface on undirected graphs
 */

#include <iostream>
#include <memory>
#include <thread>
#include <unordered_set>
#include <vector>

#include "gtest/gtest.h"
#include "third-party/libcuckoo/cuckoohash_map.hh"

#include "configuration.hpp"
#include "graph/edge_stream.hpp"
#include "library/interface.hpp"
#include "library/baseline/adjacency_list.hpp"

#if defined(HAVE_LLAMA)
#include "library/llama/llama_class.hpp"
#endif

#if defined(HAVE_STINGER)
#include "library/stinger/stinger.hpp"
#endif

#if defined(HAVE_GRAPHONE)
#include "library/graphone/graphone.hpp"
#endif

using namespace gfe::library;
using namespace std;

static std::unique_ptr<gfe::graph::WeightedEdgeStream> generate_edge_stream(uint64_t max_vector_id){
    vector<gfe::graph::WeightedEdge> edges;
    for(uint64_t i = 1; i < max_vector_id; i++){
        for(uint64_t j = i + 2; j < max_vector_id; j+=2){
            edges.push_back(gfe::graph::WeightedEdge{i, j, static_cast<double>(j * 1000 + i)});
        }
    }
    return make_unique<gfe::graph::WeightedEdgeStream>(edges);
}

static void sequential(shared_ptr<UpdateInterface> interface, bool edge_deletions = true, bool vertex_deletions = true, uint64_t num_vertices = 1024) {
    interface->on_main_init(1);
    interface->on_thread_init(0);

    // insert all edges
    unordered_set<uint64_t> vertices_contained;
    auto edge_list = generate_edge_stream(num_vertices);
    edge_list->permute();
    for(uint64_t i =0, sz = edge_list->num_edges(); i < sz; i++){
        auto edge = edge_list->get(i);
        auto p1 = vertices_contained.insert(edge.source());
        if(p1.second) interface->add_vertex(edge.source());
        auto p2 = vertices_contained.insert(edge.destination());
        if(p2.second) interface->add_vertex(edge.destination());

        // insert sometimes as <i, j> and sometimes as <j, i>
        if((edge.m_source + edge.m_destination) % 2 == 0){
            swap(edge.m_source, edge.m_destination);
        }

        interface->add_edge(edge);
    }

    interface->build();

    // check all edges have been inserted
    ASSERT_EQ(interface->num_edges(), edge_list->num_edges());
    for(uint64_t i = 1; i < edge_list->max_vertex_id(); i++){
        for(uint64_t j = i +1; j < edge_list->max_vertex_id(); j++){
            if((i + j) % 2 == 0){ // the edge should be present
                ASSERT_TRUE(interface->has_edge(i, j));
                ASSERT_TRUE(interface->has_edge(j, i)); // because it's undirected
                uint32_t expected_value = j * 1000 + i;
                ASSERT_EQ(interface->get_weight(i, j), expected_value);
                ASSERT_EQ(interface->get_weight(j, i), expected_value);
            } else {
                ASSERT_FALSE(interface->has_edge(i, j));
                ASSERT_FALSE(interface->has_edge(j, i));
            }
        }
    }

    if(edge_deletions){ // remove all edges from the graph
        edge_list->permute(gfe::configuration().seed() + 99942341);
        for(uint64_t i =0, sz = edge_list->num_edges(); i < sz; i++){
            auto edge = edge_list->get(i).edge();

            // remove sometimes as <i, j> and sometimes as <j, i>
            if((edge.m_source + edge.m_destination) % 3 == 0){
                swap(edge.m_source, edge.m_destination);
            }

            interface->remove_edge(edge);
        }

        interface->build(); // flush all deletions in LLAMA
        ASSERT_EQ(interface->num_edges(), 0);

        for(uint64_t i = 1; i < edge_list->max_vertex_id(); i++){
            for(uint64_t j = i +1; j < edge_list->max_vertex_id(); j++){
                ASSERT_FALSE(interface->has_edge(i, j));
                ASSERT_FALSE(interface->has_edge(j, i));
            }
        }

        // remove all vertices
        if(vertex_deletions){
            for(uint64_t vertex_id = 1; vertex_id <= edge_list->max_vertex_id(); vertex_id++){
                interface->remove_vertex(vertex_id);
            }
            interface->build(); // flush all changes
            ASSERT_EQ(interface->num_vertices(), 0);
        }
    }

    // done
    interface->on_thread_destroy(0);
    interface->on_main_destroy();
}

// Whether execute the tests inside #parallel() with multiple threads. Useful for GraphOne as the implementation of
// get weight is particularly slow
static bool parallel_check = false;

// Whether vertex deletions are supported by the interface. All libraries support them (more or less), but GraphOne
static bool parallel_vertex_deletions = true;

static void parallel(shared_ptr<UpdateInterface> interface, uint64_t num_vertices, int num_threads = 8, bool deletions = true) {
    assert(num_threads > 0);
    interface->on_main_init(num_threads);

    // insert all edges
    cuckoohash_map<uint64_t, bool> vertices;
    shared_ptr<gfe::graph::WeightedEdgeStream> edge_list = generate_edge_stream(num_vertices);
    edge_list->permute();

//    cout << "num edges: " << edge_list->num_edges() << endl;

    auto routine_insert_edges = [interface, edge_list, &vertices](int thread_id, uint64_t start, uint64_t length){
        interface->on_thread_init(thread_id);

        for(int64_t pos = start, end = start + length; pos < end; pos++){
            auto edge = edge_list->get(pos);

            if(vertices.insert(edge.source(), true)) interface->add_vertex(edge.source());
            if(vertices.insert(edge.destination(), true)) interface->add_vertex(edge.destination());

            // insert sometimes as <i, j> and sometimes as <j, i>
            if((edge.m_source + edge.m_destination) % 2 == 0){
                swap(edge.m_source, edge.m_destination);
            }

            // the function returns true if the edge has been inserted. Repeat the loop if it cannot insert the edge as one of
            // the vertices is still being inserted by another thread
            while( ! interface->add_edge(edge) ) { /* nop */ } ;
        }

        interface->on_thread_destroy(thread_id);
    };

    uint64_t edges_per_thread = edge_list->num_edges() / num_threads;
    uint64_t odd_threads = edge_list->num_edges() % num_threads;

    vector<thread> threads;
    uint64_t start = 0;
    for(int thread_id = 0; thread_id < num_threads; thread_id ++){
        uint64_t length = edges_per_thread + (thread_id < odd_threads);
        threads.emplace_back(routine_insert_edges, thread_id, start, length);
        start += length;
    }
    for(auto& t : threads) t.join();
    threads.clear();

    interface->on_thread_init(0);
    interface->build();
//    interface->dump();
    interface->on_thread_destroy(0);

    // check all edges have been inserted
    auto routine_check_edges = [interface, edge_list](int thread_id, int num_threads){
        interface->on_thread_destroy(thread_id);
        ASSERT_EQ(interface->num_edges(), edge_list->num_edges());
        for(uint64_t i = thread_id +1; i < edge_list->max_vertex_id(); i += num_threads){
            for(uint64_t j = i +1; j < edge_list->max_vertex_id(); j++){
                if((i + j) % 2 == 0){ // the edge should be present
                    ASSERT_TRUE(interface->has_edge(i, j));
                    ASSERT_TRUE(interface->has_edge(j, i));
                    uint32_t expected_value = j * 1000 + i;
                    ASSERT_EQ(interface->get_weight(i, j), expected_value);
                    ASSERT_EQ(interface->get_weight(j, i), expected_value);
                } else {
                    ASSERT_FALSE(interface->has_edge(i, j));
                    ASSERT_FALSE(interface->has_edge(j, i));
                }
            }
        }
        interface->on_thread_destroy(thread_id);
    };

    int num_threads_check = parallel_check ? num_threads : /* sequential */ 1;
    for(int thread_id = 0; thread_id < num_threads_check; thread_id ++){
        threads.emplace_back(routine_check_edges, thread_id, num_threads_check);
    }
    for(auto& t : threads) t.join();
    threads.clear();


    if(deletions){ // remove all edges from the graph
        edge_list->permute(gfe::configuration().seed() + 99942341);

        auto routine_remove_edges = [interface, edge_list](int thread_id, uint64_t start, uint64_t length){
            interface->on_thread_init(thread_id);

            for(int64_t pos = start, end = start + length; pos < end; pos++){
                auto edge = edge_list->get(pos).edge();

                // remove sometimes as <i, j> and sometimes as <j, i>
                if((edge.m_source + edge.m_destination) % 3 == 0){
                    swap(edge.m_source, edge.m_destination);
                }

                interface->remove_edge(edge);
            }

            interface->on_thread_destroy(thread_id);
        };
        threads.clear();
        start = 0;
        for(int thread_id = 0; thread_id < num_threads; thread_id ++){
            uint64_t length = edges_per_thread + (thread_id < odd_threads);
            threads.emplace_back(routine_remove_edges, thread_id, start, length);
            start += length;
        }
        for(auto& t : threads) t.join();

        // check all edges have been removed
        interface->on_thread_init(0);
        interface->build();
        ASSERT_EQ(interface->num_edges(), 0);
        for(uint64_t i = 1; i < edge_list->max_vertex_id(); i++){
            for(uint64_t j = i +1; j < edge_list->max_vertex_id(); j++){
                ASSERT_FALSE(interface->has_edge(i, j));
                ASSERT_FALSE(interface->has_edge(j, i));
            }
        }

        // remove all vertices from the graph
        if(parallel_vertex_deletions){
            auto routine_remove_vertices = [interface, edge_list](int thread_id, uint64_t start, uint64_t length){
                interface->on_thread_init(thread_id);

                for(int64_t pos = start, end = start + length; pos < end; pos++){
                    interface->remove_vertex(pos +1);
                }

                interface->on_thread_destroy(thread_id);
            };
            uint64_t vertices_per_thread = edge_list->max_vertex_id() / num_threads;
            odd_threads = edge_list->max_vertex_id() % num_threads;
            threads.clear();
            start = 0;
            for(int thread_id = 0; thread_id < num_threads; thread_id ++){
                uint64_t length = vertices_per_thread + (thread_id < odd_threads);
                threads.emplace_back(routine_remove_vertices, thread_id, start, length);
                start += length;
            }
            for(auto& t : threads) t.join();

            interface->build();
            ASSERT_EQ(interface->num_vertices(), 0);
        } // parallel_vertex_deletions
        interface->on_thread_destroy(0);
    }

    // done
    interface->on_main_destroy();
}

TEST(AdjacencyList, UpdatesUndirected){
    auto adjlist = make_shared<AdjacencyList>(/* directed */ false);
    sequential(adjlist);
    parallel(adjlist, 128);
    parallel(adjlist, 1024);
}

#if defined(HAVE_LLAMA)
TEST(LLAMA, UpdatesUndirected){
    auto llama = make_shared<LLAMAClass>(/* directed */ false);
    sequential(llama, /* perform deletions ? */ false);
    llama = make_shared<LLAMAClass>(/* directed */ false); // reinit as we didn't perform the deletions in `sequential'
    parallel(llama, 128);
    parallel(llama, 1024);
}
#endif

#if defined(HAVE_STINGER)
TEST(Stinger, UpdatesUndirected) {
    auto stinger = make_shared<Stinger>(/* directed */ false);
    sequential(stinger);
    parallel(stinger, 128);
    parallel(stinger, 1024);
}
#endif

#if defined(HAVE_GRAPHONE)
TEST(GraphOne, UpdatesUndirected){
    parallel_check = true; // global, check the weights in parallel
    parallel_vertex_deletions = false; // global, disable vertex deletions

    auto graphone = make_shared<GraphOne>(/* directed ? */ false, /* vtx dict ? */ true, /* blind writes ? */ false, /* ignore build = */ false, /* GAP impl ? */ false, 1ull << 24);
    sequential(graphone, true, false, 16);
    parallel(graphone, 32);

    parallel_check = false; // global, reset to the default value
    parallel_vertex_deletions = true; // global, reset to the default value
}
#endif
