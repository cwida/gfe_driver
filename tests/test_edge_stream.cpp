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

#include <iostream>
#include "common/error.hpp"
#include "common/filesystem.hpp"
#include "graph/edge.hpp"
#include "graph/edge_stream.hpp"

using namespace gfe::graph;
using namespace std;

TEST(EdgeStream, Sanity) {
    WeightedEdgeStream stream(common::filesystem::directory_executable() + "/graphs/weighted_no_comments.wel");

    ASSERT_EQ(stream.num_edges(), 10);
    ASSERT_EQ(stream.max_vertex_id(), 4);

    // 1 2 10
    ASSERT_EQ(stream[0].source(), 1);
    ASSERT_EQ(stream[0].destination(), 2);
    ASSERT_EQ(stream[0].weight(), 10);

    // 1 3 100
    ASSERT_EQ(stream[1].source(), 1);
    ASSERT_EQ(stream[1].destination(), 3);
    ASSERT_EQ(stream[1].weight(), 100);

    // 1 4 200
    ASSERT_EQ(stream[2].source(), 1);
    ASSERT_EQ(stream[2].destination(), 4);
    ASSERT_EQ(stream[2].weight(), 200);

    // 2 1 10
    ASSERT_EQ(stream[3].source(), 2);
    ASSERT_EQ(stream[3].destination(), 1);
    ASSERT_EQ(stream[3].weight(), 10);

    // 2 4 10
    ASSERT_EQ(stream[4].source(), 2);
    ASSERT_EQ(stream[4].destination(), 4);
    ASSERT_EQ(stream[4].weight(), 10);

    // 3 1 100
    ASSERT_EQ(stream[5].source(), 3);
    ASSERT_EQ(stream[5].destination(), 1);
    ASSERT_EQ(stream[5].weight(), 100);

    // 3 4 10
    ASSERT_EQ(stream[6].source(), 3);
    ASSERT_EQ(stream[6].destination(), 4);
    ASSERT_EQ(stream[6].weight(), 10);

    // 4 1 200
    ASSERT_EQ(stream[7].source(), 4);
    ASSERT_EQ(stream[7].destination(), 1);
    ASSERT_EQ(stream[7].weight(), 200);

    // 4 2 10
    ASSERT_EQ(stream[8].source(), 4);
    ASSERT_EQ(stream[8].destination(), 2);
    ASSERT_EQ(stream[8].weight(), 10);

    // 4 3 10
    ASSERT_EQ(stream[9].source(), 4);
    ASSERT_EQ(stream[9].destination(), 3);
    ASSERT_EQ(stream[9].weight(), 10);

    // permutation
    WeightedEdgeStream permuted(common::filesystem::directory_executable() + "/graphs/weighted_no_comments.wel");
    permuted.permute();
    ASSERT_EQ(stream.num_edges(), permuted.num_edges());
    for(int i = 0, end = stream.num_edges(); i < end; i++){
        bool found = false;
        for(int j = 0; j < end && !found; j++){
            found = (stream[i] == permuted[j]);
        }
        ASSERT_EQ(found, true);
    }
}


