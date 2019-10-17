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
#include "reader/graphalytics_reader.hpp"
#include "reader/metis_reader.hpp"
#include "reader/plain_reader.hpp"

using namespace std;
using namespace graph;
using namespace reader;


// validate the next read from a reader
static void validate_read(Reader& reader, uint64_t source, uint64_t dest, double weight){
    WeightedEdge edge;
    bool rc = reader.read(edge);
    ASSERT_TRUE(rc);
    ASSERT_EQ(edge.source(), source);
    ASSERT_EQ(edge.destination(), dest);
    ASSERT_EQ(edge.weight(), weight);
}

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

    // vertex 1
    // 2 10 3 100 4 200
    validate_read(reader, 1, 2, 10);
    validate_read(reader, 1, 3, 100);
    validate_read(reader, 1, 4, 200);

    // vertex 2
    // 1 10 4 10
    validate_read(reader, 2, 1, 10);
    validate_read(reader, 2, 4, 10);

    // vertex 3
    // 1 100 4 10
    validate_read(reader, 3, 1, 100);
    validate_read(reader, 3, 4, 10);

    // vertex 4
    // 1 200 2 10 3 10
    validate_read(reader, 4, 1, 200);
    validate_read(reader, 4, 2, 10);
    validate_read(reader, 4, 3, 10);

    // eof
    WeightedEdge edge;
    ASSERT_FALSE( reader.read(edge) );
}

TEST(Metis, WithComments) {
    // same graph as metis without comments, the parser should be able to skip all commented lines successfully

    MetisReader reader(common::filesystem::directory_executable() + "/graphs/weighted_with_comments.metis");

    // vertex 1
    // 2 10 3 100 4 200
    validate_read(reader, 1, 2, 10);
    validate_read(reader, 1, 3, 100);
    validate_read(reader, 1, 4, 200);

    // vertex 2
    // 1 10 4 10
    validate_read(reader, 2, 1, 10);
    validate_read(reader, 2, 4, 10);

    // vertex 3
    // 1 100 4 10
    validate_read(reader, 3, 1, 100);
    validate_read(reader, 3, 4, 10);

    // vertex 4
    // 1 200 2 10 3 10
    validate_read(reader, 4, 1, 200);
    validate_read(reader, 4, 2, 10);
    validate_read(reader, 4, 3, 10);

    // eof
    WeightedEdge edge;
    ASSERT_FALSE( reader.read(edge) );
}

TEST(Dimacs9, WithoutComments) {
    Dimacs9Reader reader(common::filesystem::directory_executable() + "/graphs/weighted_no_comments.dimacs9");

    // vertex 1
    // 2 10 3 100 4 200
    validate_read(reader, 1, 2, 10);
    validate_read(reader, 1, 3, 100);
    validate_read(reader, 1, 4, 200);

    // vertex 2
    // 1 10 4 10
    validate_read(reader, 2, 1, 10);
    validate_read(reader, 2, 4, 10);

    // vertex 3
    // 1 100 4 10
    validate_read(reader, 3, 1, 100);
    validate_read(reader, 3, 4, 10);

    // vertex 4
    // 1 200 2 10 3 10
    validate_read(reader, 4, 1, 200);
    validate_read(reader, 4, 2, 10);
    validate_read(reader, 4, 3, 10);

    // eof
    WeightedEdge edge;
    ASSERT_FALSE( reader.read(edge) );
}

TEST(Dimacs9, WithComments) {
    Dimacs9Reader reader(common::filesystem::directory_executable() + "/graphs/weighted_with_comments.dimacs9");

    // vertex 1
    // 2 10 3 100 4 200
    validate_read(reader, 1, 2, 10);
    validate_read(reader, 1, 3, 100);
    validate_read(reader, 1, 4, 200);

    // vertex 2
    // 1 10 4 10
    validate_read(reader, 2, 1, 10);
    validate_read(reader, 2, 4, 10);

    // vertex 3
    // 1 100 4 10
    validate_read(reader, 3, 1, 100);
    validate_read(reader, 3, 4, 10);

    // vertex 4
    // 1 200 2 10 3 10
    validate_read(reader, 4, 1, 200);
    validate_read(reader, 4, 2, 10);
    validate_read(reader, 4, 3, 10);

    // eof
    WeightedEdge edge;
    ASSERT_FALSE( reader.read(edge) );
}

TEST(Dimacs9, Rome99) {
    Dimacs9Reader reader(common::filesystem::directory_executable() + "/graphs/rome99.gr");

    // check explicitly the first two edges
    // a 2958 2959 535
    // a 1494 1573 77
    validate_read(reader, 2958, 2959, 535);
    validate_read(reader, 1494, 1573, 77);

    // there are still 8868 in the file, we'll explicitly check only the last two
    for(int i = 0; i < 8866; i++){
        WeightedEdge edge;
        bool rc = reader.read(edge);
        ASSERT_EQ(rc, true);
    }

    // the last two edges in the graph
    // a 1545 1540 119
    // a 1374 3353 211
    validate_read(reader, 1545, 1540, 119);
    validate_read(reader, 1374, 3353, 211);

    // eof
    WeightedEdge edge;
    ASSERT_FALSE( reader.read(edge) );
}


TEST(Graphalytics, ExampleDirected) {
    GraphalyticsReader reader(common::filesystem::directory_executable() + "/graphs/ldbc_graphalytics/example-directed.properties");

    // Read the edges
    validate_read(reader, 1, 3, 0.5);
    validate_read(reader, 1, 5, 0.3);
    validate_read(reader, 2, 4, 0.1);
    validate_read(reader, 2, 5, 0.3);

    // Restart
    reader.reset();
    validate_read(reader, 1, 3, 0.5);
    validate_read(reader, 1, 5, 0.3);
    validate_read(reader, 2, 4, 0.1);
    validate_read(reader, 2, 5, 0.3);
    validate_read(reader, 2, 10, 0.12);
    validate_read(reader, 3, 1, 0.53);
    validate_read(reader, 3, 5, 0.62);
    validate_read(reader, 3, 8, 0.21);
    validate_read(reader, 3, 10, 0.52);
    validate_read(reader, 5, 3, 0.69);
    validate_read(reader, 5, 4, 0.53);
    validate_read(reader, 5, 8, 0.1);
    validate_read(reader, 6, 3, 0.23);
    validate_read(reader, 6, 4, 0.39);
    validate_read(reader, 7, 4, 0.83);
    validate_read(reader, 8, 1, 0.39);
    validate_read(reader, 9, 4, 0.69);

    // Done
    WeightedEdge edge;
    ASSERT_FALSE(reader.read(edge));

    // Read the vertices
    uint64_t vertex { 0 };
    ASSERT_TRUE( reader.read_vertex(vertex) );
    ASSERT_EQ(vertex, 1);
    ASSERT_TRUE( reader.read_vertex(vertex) );
    ASSERT_EQ(vertex, 2);
    ASSERT_TRUE( reader.read_vertex(vertex) );
    ASSERT_EQ(vertex, 3);
    ASSERT_TRUE( reader.read_vertex(vertex) );
    ASSERT_EQ(vertex, 4);
    ASSERT_TRUE( reader.read_vertex(vertex) );
    ASSERT_EQ(vertex, 5);
    ASSERT_TRUE( reader.read_vertex(vertex) );
    ASSERT_EQ(vertex, 6);
    ASSERT_TRUE( reader.read_vertex(vertex) );
    ASSERT_EQ(vertex, 7);
    ASSERT_TRUE( reader.read_vertex(vertex) );
    ASSERT_EQ(vertex, 8);
    ASSERT_TRUE( reader.read_vertex(vertex) );
    ASSERT_EQ(vertex, 9);
    ASSERT_TRUE( reader.read_vertex(vertex) );
    ASSERT_EQ(vertex, 10);
}

TEST(Graphalytics, ExampleUndirected) {
    GraphalyticsReader reader(common::filesystem::directory_executable() + "/graphs/ldbc_graphalytics/example-undirected.properties");

    // Read the edges
    reader.set_emit_directed_edges(true);
    validate_read(reader, 2, 3, 0.9);
    validate_read(reader, 3, 2, 0.9);
    validate_read(reader, 2, 4, 0.69);
    validate_read(reader, 4, 2, 0.69);

    // Restart
    reader.reset();
    reader.set_emit_directed_edges(false);
    validate_read(reader, 2, 3, 0.9);
    validate_read(reader, 2, 4, 0.69);
    validate_read(reader, 3, 4, 0.13);
    validate_read(reader, 3, 5, 0.5);
    validate_read(reader, 3, 8, 0.32);
    validate_read(reader, 5, 6, 0.63);
    validate_read(reader, 5, 8, 0.12);
    validate_read(reader, 6, 7, 0.53);
    validate_read(reader, 6, 8, 0.64);
    validate_read(reader, 6, 9, 0.23);
    validate_read(reader, 6, 10, 0.63);
    validate_read(reader, 7, 9, 0.36);
}
