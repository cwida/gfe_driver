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
#include "reader/dimacs9_reader.hpp"
#include "reader/metis_reader.hpp"
#include "reader/plain_reader.hpp"

using namespace std;
using namespace reader;

TEST(PlainWeighted, WithoutComments) {
    graph::WeightedEdgeStream stream{  common::filesystem::directory_executable() + "/graphs/weighted_no_comments.wel" };
    stream.permute();

}

TEST(PlainWeighted, WithComments) {
    // same graph as the one with no comments, the parser should be able to skip all commented lines successfully
    PlainReader reader(common::filesystem::directory_executable() + "/graphs/weighted_with_comments.wel", /* weighted */ true);

    WeightedEdge edge;
    bool rc { false };

    // vertex 1
    // 2 10 3 100 4 200
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 1);
    ASSERT_EQ(edge.destination(), 2);
    ASSERT_EQ(edge.weight(), 10);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 1);
    ASSERT_EQ(edge.destination(), 3);
    ASSERT_EQ(edge.weight(), 100);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 1);
    ASSERT_EQ(edge.destination(), 4);
    ASSERT_EQ(edge.weight(), 200);

    // vertex 2
    // 1 10 4 10
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 2);
    ASSERT_EQ(edge.destination(), 1);
    ASSERT_EQ(edge.weight(), 10);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 2);
    ASSERT_EQ(edge.destination(), 4);
    ASSERT_EQ(edge.weight(), 10);

    // vertex 3
    // 1 100 4 10
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 3);
    ASSERT_EQ(edge.destination(), 1);
    ASSERT_EQ(edge.weight(), 100);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 3);
    ASSERT_EQ(edge.destination(), 4);
    ASSERT_EQ(edge.weight(), 10);

    // vertex 4
    // 1 200 2 10 3 10
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 4);
    ASSERT_EQ(edge.destination(), 1);
    ASSERT_EQ(edge.weight(), 200);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 4);
    ASSERT_EQ(edge.destination(), 2);
    ASSERT_EQ(edge.weight(), 10);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 4);
    ASSERT_EQ(edge.destination(), 3);
    ASSERT_EQ(edge.weight(), 10);

    // eof
    rc = reader.read(edge);
    ASSERT_EQ(rc, false);
}

TEST(Plain, ConstantWeight) {
    PlainReader reader(common::filesystem::directory_executable() + "/graphs/non_weighted.el", /* weighted */ false, /* max weight */ 1);

    WeightedEdge edge;
    bool rc { false };

    // vertex 1
    // 2 10 3 100 4 200
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 1);
    ASSERT_EQ(edge.destination(), 2);
    ASSERT_EQ(edge.weight(), 1);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 1);
    ASSERT_EQ(edge.destination(), 3);
    ASSERT_EQ(edge.weight(), 1);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 1);
    ASSERT_EQ(edge.destination(), 4);
    ASSERT_EQ(edge.weight(), 1);

    // vertex 2
    // 1 10 4 10
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 2);
    ASSERT_EQ(edge.destination(), 1);
    ASSERT_EQ(edge.weight(), 1);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 2);
    ASSERT_EQ(edge.destination(), 4);
    ASSERT_EQ(edge.weight(), 1);

    // vertex 3
    // 1 100 4 10
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 3);
    ASSERT_EQ(edge.destination(), 1);
    ASSERT_EQ(edge.weight(), 1);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 3);
    ASSERT_EQ(edge.destination(), 4);
    ASSERT_EQ(edge.weight(), 1);

    // vertex 4
    // 1 200 2 10 3 10
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 4);
    ASSERT_EQ(edge.destination(), 1);
    ASSERT_EQ(edge.weight(), 1);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 4);
    ASSERT_EQ(edge.destination(), 2);
    ASSERT_EQ(edge.weight(), 1);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 4);
    ASSERT_EQ(edge.destination(), 3);
    ASSERT_EQ(edge.weight(), 1);

    // eof
    rc = reader.read(edge);
    ASSERT_EQ(rc, false);
}

TEST(Plain, RandomWeight) {
    PlainReader reader(common::filesystem::directory_executable() + "/graphs/non_weighted.el", /* weighted */ false, /* max weight */ 5);

    WeightedEdge edge;
    bool rc { false };

    // vertex 1
    // 2 10 3 100 4 200
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 1);
    ASSERT_EQ(edge.destination(), 2);
    ASSERT_GE(edge.weight(), 1); ASSERT_LE(edge.weight(), 5); // weight in [1, 5]
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 1);
    ASSERT_EQ(edge.destination(), 3);
    ASSERT_GE(edge.weight(), 1); ASSERT_LE(edge.weight(), 5); // weight in [1, 5]
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 1);
    ASSERT_EQ(edge.destination(), 4);
    ASSERT_GE(edge.weight(), 1); ASSERT_LE(edge.weight(), 5); // weight in [1, 5]

    // vertex 2
    // 1 10 4 10
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 2);
    ASSERT_EQ(edge.destination(), 1);
    ASSERT_GE(edge.weight(), 1); ASSERT_LE(edge.weight(), 5); // weight in [1, 5]
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 2);
    ASSERT_EQ(edge.destination(), 4);
    ASSERT_GE(edge.weight(), 1); ASSERT_LE(edge.weight(), 5); // weight in [1, 5]

    // vertex 3
    // 1 100 4 10
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 3);
    ASSERT_EQ(edge.destination(), 1);
    ASSERT_GE(edge.weight(), 1); ASSERT_LE(edge.weight(), 5); // weight in [1, 5]
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 3);
    ASSERT_EQ(edge.destination(), 4);
    ASSERT_GE(edge.weight(), 1); ASSERT_LE(edge.weight(), 5); // weight in [1, 5]

    // vertex 4
    // 1 200 2 10 3 10
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 4);
    ASSERT_EQ(edge.destination(), 1);
    ASSERT_GE(edge.weight(), 1); ASSERT_LE(edge.weight(), 5); // weight in [1, 5]
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 4);
    ASSERT_EQ(edge.destination(), 2);
    ASSERT_GE(edge.weight(), 1); ASSERT_LE(edge.weight(), 5); // weight in [1, 5]
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 4);
    ASSERT_EQ(edge.destination(), 3);
    ASSERT_GE(edge.weight(), 1); ASSERT_LE(edge.weight(), 5); // weight in [1, 5]

    // eof
    rc = reader.read(edge);
    ASSERT_EQ(rc, false);
}

TEST(Metis, WithoutComments) {
    MetisReader reader(common::filesystem::directory_executable() + "/graphs/weighted_no_comments.metis");

    WeightedEdge edge;
    bool rc { false };

    // vertex 1
    // 2 10 3 100 4 200
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 1);
    ASSERT_EQ(edge.destination(), 2);
    ASSERT_EQ(edge.weight(), 10);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 1);
    ASSERT_EQ(edge.destination(), 3);
    ASSERT_EQ(edge.weight(), 100);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 1);
    ASSERT_EQ(edge.destination(), 4);
    ASSERT_EQ(edge.weight(), 200);

    // vertex 2
    // 1 10 4 10
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 2);
    ASSERT_EQ(edge.destination(), 1);
    ASSERT_EQ(edge.weight(), 10);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 2);
    ASSERT_EQ(edge.destination(), 4);
    ASSERT_EQ(edge.weight(), 10);

    // vertex 3
    // 1 100 4 10
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 3);
    ASSERT_EQ(edge.destination(), 1);
    ASSERT_EQ(edge.weight(), 100);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 3);
    ASSERT_EQ(edge.destination(), 4);
    ASSERT_EQ(edge.weight(), 10);

    // vertex 4
    // 1 200 2 10 3 10
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 4);
    ASSERT_EQ(edge.destination(), 1);
    ASSERT_EQ(edge.weight(), 200);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 4);
    ASSERT_EQ(edge.destination(), 2);
    ASSERT_EQ(edge.weight(), 10);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 4);
    ASSERT_EQ(edge.destination(), 3);
    ASSERT_EQ(edge.weight(), 10);

    // eof
    rc = reader.read(edge);
    ASSERT_EQ(rc, false);
}

TEST(Metis, WithComments) {
    // same graph as metis without comments, the parser should be able to skip all commented lines successfully

    MetisReader reader(common::filesystem::directory_executable() + "/graphs/weighted_with_comments.metis");

    WeightedEdge edge;
    bool rc { false };

    // vertex 1
    // 2 10 3 100 4 200
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 1);
    ASSERT_EQ(edge.destination(), 2);
    ASSERT_EQ(edge.weight(), 10);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 1);
    ASSERT_EQ(edge.destination(), 3);
    ASSERT_EQ(edge.weight(), 100);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 1);
    ASSERT_EQ(edge.destination(), 4);
    ASSERT_EQ(edge.weight(), 200);

    // vertex 2
    // 1 10 4 10
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 2);
    ASSERT_EQ(edge.destination(), 1);
    ASSERT_EQ(edge.weight(), 10);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 2);
    ASSERT_EQ(edge.destination(), 4);
    ASSERT_EQ(edge.weight(), 10);

    // vertex 3
    // 1 100 4 10
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 3);
    ASSERT_EQ(edge.destination(), 1);
    ASSERT_EQ(edge.weight(), 100);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 3);
    ASSERT_EQ(edge.destination(), 4);
    ASSERT_EQ(edge.weight(), 10);

    // vertex 4
    // 1 200 2 10 3 10
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 4);
    ASSERT_EQ(edge.destination(), 1);
    ASSERT_EQ(edge.weight(), 200);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 4);
    ASSERT_EQ(edge.destination(), 2);
    ASSERT_EQ(edge.weight(), 10);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 4);
    ASSERT_EQ(edge.destination(), 3);
    ASSERT_EQ(edge.weight(), 10);

    // eof
    rc = reader.read(edge);
    ASSERT_EQ(rc, false);
}

TEST(Dimacs9, WithoutComments) {
    Dimacs9Reader reader(common::filesystem::directory_executable() + "/graphs/weighted_no_comments.dimacs9");

    WeightedEdge edge;
    bool rc { false };

    // vertex 1
    // 2 10 3 100 4 200
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 1);
    ASSERT_EQ(edge.destination(), 2);
    ASSERT_EQ(edge.weight(), 10);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 1);
    ASSERT_EQ(edge.destination(), 3);
    ASSERT_EQ(edge.weight(), 100);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 1);
    ASSERT_EQ(edge.destination(), 4);
    ASSERT_EQ(edge.weight(), 200);

    // vertex 2
    // 1 10 4 10
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 2);
    ASSERT_EQ(edge.destination(), 1);
    ASSERT_EQ(edge.weight(), 10);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 2);
    ASSERT_EQ(edge.destination(), 4);
    ASSERT_EQ(edge.weight(), 10);

    // vertex 3
    // 1 100 4 10
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 3);
    ASSERT_EQ(edge.destination(), 1);
    ASSERT_EQ(edge.weight(), 100);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 3);
    ASSERT_EQ(edge.destination(), 4);
    ASSERT_EQ(edge.weight(), 10);

    // vertex 4
    // 1 200 2 10 3 10
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 4);
    ASSERT_EQ(edge.destination(), 1);
    ASSERT_EQ(edge.weight(), 200);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 4);
    ASSERT_EQ(edge.destination(), 2);
    ASSERT_EQ(edge.weight(), 10);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 4);
    ASSERT_EQ(edge.destination(), 3);
    ASSERT_EQ(edge.weight(), 10);

    // eof
    rc = reader.read(edge);
    ASSERT_EQ(rc, false);
}

TEST(Dimacs9, WithComments) {
    Dimacs9Reader reader(common::filesystem::directory_executable() + "/graphs/weighted_with_comments.dimacs9");

    WeightedEdge edge;
    bool rc { false };

    // vertex 1
    // 2 10 3 100 4 200
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 1);
    ASSERT_EQ(edge.destination(), 2);
    ASSERT_EQ(edge.weight(), 10);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 1);
    ASSERT_EQ(edge.destination(), 3);
    ASSERT_EQ(edge.weight(), 100);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 1);
    ASSERT_EQ(edge.destination(), 4);
    ASSERT_EQ(edge.weight(), 200);

    // vertex 2
    // 1 10 4 10
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 2);
    ASSERT_EQ(edge.destination(), 1);
    ASSERT_EQ(edge.weight(), 10);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 2);
    ASSERT_EQ(edge.destination(), 4);
    ASSERT_EQ(edge.weight(), 10);

    // vertex 3
    // 1 100 4 10
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 3);
    ASSERT_EQ(edge.destination(), 1);
    ASSERT_EQ(edge.weight(), 100);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 3);
    ASSERT_EQ(edge.destination(), 4);
    ASSERT_EQ(edge.weight(), 10);

    // vertex 4
    // 1 200 2 10 3 10
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 4);
    ASSERT_EQ(edge.destination(), 1);
    ASSERT_EQ(edge.weight(), 200);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 4);
    ASSERT_EQ(edge.destination(), 2);
    ASSERT_EQ(edge.weight(), 10);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 4);
    ASSERT_EQ(edge.destination(), 3);
    ASSERT_EQ(edge.weight(), 10);

    // eof
    rc = reader.read(edge);
    ASSERT_EQ(rc, false);
}

TEST(Dimacs9, Rome99) {
    Dimacs9Reader reader(common::filesystem::directory_executable() + "/graphs/rome99.gr");

    WeightedEdge edge;
    bool rc { false };

    // check explicitly the first two edges
    // a 2958 2959 535
    // a 1494 1573 77
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 2958);
    ASSERT_EQ(edge.destination(), 2959);
    ASSERT_EQ(edge.weight(), 535);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 1494);
    ASSERT_EQ(edge.destination(), 1573);
    ASSERT_EQ(edge.weight(), 77);

    // there are still 8868 in the file, we'll explicitly check only the last two
    for(int i = 0; i < 8866; i++){
        rc = reader.read(edge);
        ASSERT_EQ(rc, true);
    }

    // the last two edges in the graph
    // a 1545 1540 119
    // a 1374 3353 211
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 1545);
    ASSERT_EQ(edge.destination(), 1540);
    ASSERT_EQ(edge.weight(), 119);
    rc = reader.read(edge);
    ASSERT_EQ(rc, true);
    ASSERT_EQ(edge.source(), 1374);
    ASSERT_EQ(edge.destination(), 3353);
    ASSERT_EQ(edge.weight(), 211);

    // eof
    rc = reader.read(edge);
    ASSERT_EQ(rc, false);
}
