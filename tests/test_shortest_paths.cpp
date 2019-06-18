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

#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>
#include <unordered_set>

#include "common/system.hpp"
#include "gtest/gtest.h"
#include "graph/edge_stream.hpp"
#include "library/interface.hpp"
#include "configuration.hpp"

#if defined(HAVE_STINGER)
#include "library/stinger/stinger.hpp"
#endif

using namespace library;
using namespace std;

static void load(shared_ptr<UpdateInterface> interface){
    ASSERT_EQ(interface->num_edges(), 0); // the interface is not empty!

    graph::WeightedEdgeStream stream{  common::filesystem::directory_executable() + "/graphs/shortest_paths.wel" };
    stream.permute();
    unordered_set<uint64_t> vertices;

    for(size_t i = 0, sz = stream.num_edges(); i < sz; i++){
        auto edge = stream.get(i);
        if(vertices.insert(edge.source()).second) interface->add_vertex(edge.source());
        if(vertices.insert(edge.destination()).second) interface->add_vertex(edge.destination());
        interface->add_edge(edge);
    }

    ASSERT_EQ(interface->num_edges(), stream.num_edges());
}

static void validate_bfs(shared_ptr<ShortestPathInterface> interface){
    // BFS 1 : 1
    ASSERT_EQ(interface->bfs_one(1, 1), 0 );
    ASSERT_EQ(interface->bfs_one(1, 10), 1);
    ASSERT_EQ(interface->bfs_one(1, 11), 1);
    ASSERT_EQ(interface->bfs_one(1, 12), 1);
    ASSERT_EQ(interface->bfs_one(1, 20), 2);
    ASSERT_EQ(interface->bfs_one(1, 21), 2);
    ASSERT_EQ(interface->bfs_one(1, 22), 2);
    ASSERT_EQ(interface->bfs_one(1, 23), 2);
    ASSERT_EQ(interface->bfs_one(1, 24), 2);
    ASSERT_EQ(interface->bfs_one(1, 30), 3);
    ASSERT_EQ(interface->bfs_one(1, 41), 4);
    ASSERT_EQ(interface->bfs_one(1, 42), 4);
    ASSERT_EQ(interface->bfs_one(1, 35), -1); // unreachable
    ASSERT_EQ(interface->bfs_one(1, 45), -1); // unreachable

    std::vector<ShortestPathInterface::Distance> path;
    interface->bfs_one(1, 42, &path);
    ASSERT_EQ(path.size(), 3);
    ASSERT_EQ(path.at(0).m_distance, 1);
    ASSERT_EQ(path.at(1).m_distance, 2);
    ASSERT_EQ(path.at(2).m_distance, 3);
    ASSERT_TRUE(path.at(0).m_vertex == 11 || path.at(0).m_vertex == 12 || path.at(0).m_vertex == 13);
    ASSERT_TRUE(path.at(1).m_vertex >= 20 && path.at(1).m_vertex <= 24);
    ASSERT_TRUE(path.at(2).m_vertex == 30);

    // BFS 1 : all
    std::vector<ShortestPathInterface::Distance> distances;
    interface->bfs_all(1, &distances);
    ASSERT_GE(distances.size(), 11); // it may also have reported the nodes unreacheable
    sort(begin(distances), end(distances), [](const auto& d1, const auto& d2){
        return (d1.m_distance < d2.m_distance) || (d1.m_distance == d2.m_distance && d1.m_vertex < d2.m_vertex);
    });

//    for(size_t i = 0; i < distances.size(); i++){
//        cout << "[" << i << "] " << distances[i] << "\n";
//    }

    // ignore the distance to the source
    if(distances.at(0).m_vertex == 1){
        ASSERT_EQ(distances.at(0).m_distance, 0);
        ASSERT_GE(distances.size(), 12);
        distances.erase(begin(distances));
    }
    ASSERT_EQ(distances[0].m_vertex, 10);
    ASSERT_EQ(distances[0].m_distance, 1);
    ASSERT_EQ(distances[1].m_vertex, 11);
    ASSERT_EQ(distances[1].m_distance, 1);
    ASSERT_EQ(distances[2].m_vertex, 12);
    ASSERT_EQ(distances[2].m_distance, 1);
    ASSERT_EQ(distances[3].m_vertex, 20);
    ASSERT_EQ(distances[3].m_distance, 2);
    ASSERT_EQ(distances[4].m_vertex, 21);
    ASSERT_EQ(distances[4].m_distance, 2);
    ASSERT_EQ(distances[5].m_vertex, 22);
    ASSERT_EQ(distances[5].m_distance, 2);
    ASSERT_EQ(distances[6].m_vertex, 23);
    ASSERT_EQ(distances[6].m_distance, 2);
    ASSERT_EQ(distances[7].m_vertex, 24);
    ASSERT_EQ(distances[7].m_distance, 2);
    ASSERT_EQ(distances[8].m_vertex, 30);
    ASSERT_EQ(distances[8].m_distance, 3);
    ASSERT_EQ(distances[9].m_vertex, 41);
    ASSERT_EQ(distances[9].m_distance, 4);
    ASSERT_EQ(distances[10].m_vertex, 42);
    ASSERT_EQ(distances[10].m_distance, 4);

    // in case also the unreacheable vertices have been reported, expect an infinite distance
    for(size_t i = 11; i < distances.size(); i++){
        ASSERT_GE(distances[i].m_distance, std::numeric_limits<int64_t>::max());
    }
}

static void validate_spw(shared_ptr<ShortestPathInterface> interface){
    // Shortest paths 1 : 1
    ASSERT_EQ(interface->spw_one(1, 1), 0 );
    ASSERT_EQ(interface->spw_one(1, 10), 5);
    ASSERT_EQ(interface->spw_one(1, 11), 10);
    ASSERT_EQ(interface->spw_one(1, 12), 15);
    ASSERT_EQ(interface->spw_one(1, 20), 95);
    ASSERT_EQ(interface->spw_one(1, 21), 80);
    ASSERT_EQ(interface->spw_one(1, 22), 70);
    ASSERT_EQ(interface->spw_one(1, 23), 55);
    ASSERT_EQ(interface->spw_one(1, 24), 40);
    ASSERT_EQ(interface->spw_one(1, 30), 45);
    ASSERT_EQ(interface->spw_one(1, 41), 60);
    ASSERT_EQ(interface->spw_one(1, 42), 50);
    ASSERT_EQ(interface->spw_one(1, 35), -1); // unreachable
    ASSERT_EQ(interface->spw_one(1, 45), -1); // unreachable

    std::vector<ShortestPathInterface::Distance> path;
    interface->spw_one(1, 41, &path);
    ASSERT_EQ(path.size(), 4);
    ASSERT_EQ(path[0].m_vertex, 12);
    ASSERT_EQ(path[0].m_distance, 15);
    ASSERT_EQ(path[1].m_vertex, 24);
    ASSERT_EQ(path[1].m_distance, 40);
    ASSERT_EQ(path[2].m_vertex, 30);
    ASSERT_EQ(path[2].m_distance, 45);
    ASSERT_EQ(path[3].m_vertex, 42);
    ASSERT_EQ(path[3].m_distance, 50);

    // Shortest paths 1 : all
    std::vector<ShortestPathInterface::Distance> distances;
    interface->spw_all(1, &distances);
    ASSERT_GE(distances.size(), 11); // it may also have reported the nodes unreacheable
    sort(begin(distances), end(distances), [](const auto& d1, const auto& d2){
        return (d1.m_distance < d2.m_distance) || (d1.m_distance == d2.m_distance && d1.m_vertex < d2.m_vertex);
    });

//    for(size_t i = 0; i < distances.size(); i++){
//        cout << "[" << i << "] " << distances[i] << "\n";
//    }

    // ignore the distance to the source
    if(distances.at(0).m_vertex == 1){
        ASSERT_EQ(distances.at(0).m_distance, 0);
        ASSERT_GE(distances.size(), 12);
        distances.erase(begin(distances));
    }
    ASSERT_EQ(distances[0].m_vertex, 10);
    ASSERT_EQ(distances[0].m_distance, 5);
    ASSERT_EQ(distances[1].m_vertex, 11);
    ASSERT_EQ(distances[1].m_distance, 10);
    ASSERT_EQ(distances[2].m_vertex, 12);
    ASSERT_EQ(distances[2].m_distance, 15);
    ASSERT_EQ(distances[3].m_vertex, 24);
    ASSERT_EQ(distances[3].m_distance, 40);
    ASSERT_EQ(distances[4].m_vertex, 30);
    ASSERT_EQ(distances[4].m_distance, 45);
    ASSERT_EQ(distances[5].m_vertex, 42);
    ASSERT_EQ(distances[5].m_distance, 50);
    ASSERT_EQ(distances[6].m_vertex, 23);
    ASSERT_EQ(distances[6].m_distance, 55);
    ASSERT_EQ(distances[7].m_vertex, 41);
    ASSERT_EQ(distances[7].m_distance, 60);
    ASSERT_EQ(distances[8].m_vertex, 22);
    ASSERT_EQ(distances[8].m_distance, 70);
    ASSERT_EQ(distances[9].m_vertex, 21);
    ASSERT_EQ(distances[9].m_distance, 80);
    ASSERT_EQ(distances[10].m_vertex, 20);
    ASSERT_EQ(distances[10].m_distance, 95);

    // in case also the unreacheable vertices have been reported, expect an infinite distance
    for(size_t i = 11; i < distances.size(); i++){
        ASSERT_GE(distances[i].m_distance, std::numeric_limits<int64_t>::max());
    }
}

#if defined(HAVE_STINGER)
TEST(Stinger, ShortestPaths){
    auto stinger = make_shared<Stinger>();
    stinger->on_main_init(1);
    stinger->on_thread_init(0);

    load(stinger);

//    stinger->dump();

    validate_bfs(stinger);
    validate_spw(stinger);

    stinger->on_thread_destroy(0);
    stinger->on_main_destroy();
}
#endif

TEST(Dummy, ShortestPaths){
    // Placeholder in case there are no libraries installed
}
