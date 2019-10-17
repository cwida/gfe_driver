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

#include "common/filesystem.hpp"
#include "graph/edge.hpp"
#include "reader/graphlog_reader.hpp"

using namespace std;
using namespace graph;
using namespace reader;
using namespace reader::graphlog;

const std::string path_graph = common::filesystem::directory_executable() + "/graphs/ldbc_graphalytics/example-undirected.graphlog" ;

TEST(Graphlog, Properties){
    Properties properties = parse_properties(path_graph);
    ASSERT_EQ(properties["aging_coeff"], "10");
    ASSERT_EQ(properties["ef_edges"], "1");
    ASSERT_EQ(properties["ef_vertices"], "1.2");

    // check the following properties are present
    ASSERT_TRUE(properties.count("internal.edges.begin"));
    ASSERT_TRUE(properties.count("internal.edges.block_size"));
    ASSERT_TRUE(properties.count("internal.edges.cardinality"));
    ASSERT_TRUE(properties.count("internal.edges.final"));
    ASSERT_TRUE(properties.count("internal.edges.num_blocks"));
    ASSERT_TRUE(properties.count("internal.vertices.cardinality"));
    ASSERT_TRUE(properties.count("internal.vertices.final.begin"));
    ASSERT_TRUE(properties.count("internal.vertices.final.cardinality"));
    ASSERT_TRUE(properties.count("internal.vertices.temporary.begin"));
    ASSERT_TRUE(properties.count("internal.vertices.temporary.cardinality"));

    // we created this graphlog ad hoc, we already know how many vertices and edges contain
    ASSERT_EQ(properties["internal.edges.cardinality"], "120");
    ASSERT_EQ(properties["internal.edges.final"], "12");
    ASSERT_EQ(properties["internal.vertices.cardinality"], "11");
    ASSERT_EQ(properties["internal.vertices.final.cardinality"], "9");
    ASSERT_EQ(properties["internal.vertices.temporary.cardinality"], "2");
}


TEST(Graphlog, VertexLoader){
    const uint64_t vertices_sz = 4;
    unique_ptr<uint64_t[]> ptr_vertices { new uint64_t[vertices_sz] };
    uint64_t* vertices = ptr_vertices.get();

    fstream handle(path_graph, ios_base::in | ios_base::binary);
    Properties properties = parse_properties(handle);

    for(uint64_t i =0; i < 2; i++){ // repeat the test twice
        // Read the final vertices from the log
        graphlog::set_marker(properties, handle, Section::VTX_FINAL);

        VertexLoader loader { handle };
        uint64_t num_vertices = loader.load(vertices, vertices_sz);
        ASSERT_EQ(num_vertices, 4ull);
        ASSERT_EQ(vertices[0], 2ull);
        ASSERT_EQ(vertices[1], 3ull);
        ASSERT_EQ(vertices[2], 4ull);
        ASSERT_EQ(vertices[3], 5ull);
        num_vertices = loader.load(vertices, vertices_sz);
        ASSERT_EQ(num_vertices, 4ull);
        ASSERT_EQ(vertices[0], 6ull);
        ASSERT_EQ(vertices[1], 7ull);
        ASSERT_EQ(vertices[2], 8ull);
        ASSERT_EQ(vertices[3], 9ull);
        num_vertices = loader.load(vertices, vertices_sz);
        ASSERT_EQ(num_vertices, 1ull);
        ASSERT_EQ(vertices[0], 10ull);

        ASSERT_EQ( loader.load(vertices, vertices_sz), 0ull ); // depleted
        ASSERT_EQ( loader.load(vertices, vertices_sz), 0ull ); // depleted
    }

    for(uint64_t i =0; i < 2; i++) { // repeat the test twice
        // Read the temporary vertices from the log
        graphlog::set_marker(properties, handle, Section::VTX_TEMP);
        VertexLoader loader { handle };

        uint64_t num_vertices = loader.load(vertices, vertices_sz);
        ASSERT_EQ(num_vertices, 2ull);
        ASSERT_EQ(vertices[0], 1ull);
        ASSERT_EQ(vertices[1], 11ull);

        ASSERT_EQ( loader.load(vertices, vertices_sz), 0ull ); // depleted
        ASSERT_EQ( loader.load(vertices, vertices_sz), 0ull ); // depleted
    }

    handle.close();
}


TEST(Graphlog, VertexReader){
    fstream handle(path_graph, ios_base::in | ios_base::binary);
    Properties properties = parse_properties(handle);

    for(uint64_t i = 0; i < 2; i++){ // repeat the test twice
        // Read the final vertices from the log
        graphlog::set_marker(properties, handle, Section::VTX_FINAL);

        VertexReader reader { handle };
        uint64_t vertex_id = 0;

        for(uint64_t expected_vertex_id = 2; expected_vertex_id <= 10; expected_vertex_id ++){
            ASSERT_TRUE( reader.read_vertex(vertex_id) );
            ASSERT_EQ( vertex_id, expected_vertex_id );
        }

        ASSERT_FALSE( reader.read_vertex(vertex_id) ); // depleted
        ASSERT_FALSE( reader.read_vertex(vertex_id) ); // depleted
    }

    for(uint64_t i = 0; i < 2; i++){ // repeat the test twice
        // Read the final vertices from the log
        graphlog::set_marker(properties, handle, Section::VTX_TEMP);

        VertexReader reader { handle };
        uint64_t vertex_id = 0;

        ASSERT_TRUE( reader.read_vertex(vertex_id) );
        ASSERT_EQ( vertex_id, 1ull );
        ASSERT_TRUE( reader.read_vertex(vertex_id) );
        ASSERT_EQ( vertex_id, 11ul );

        ASSERT_FALSE( reader.read_vertex(vertex_id) ); // depleted
        ASSERT_FALSE( reader.read_vertex(vertex_id) ); // depleted
    }

    handle.close();
}

TEST(Graphlog, EdgeLoader){
    fstream handle(path_graph, ios_base::in | ios_base::binary);
    Properties properties = parse_properties(handle);
    graphlog::set_marker(properties, handle, Section::EDGES);

    const uint64_t array_sz = stoull(properties["internal.edges.block_size"]) / sizeof(uint64_t);
    std::unique_ptr<uint64_t[]> ptr_array { new uint64_t[array_sz] };
    uint64_t* array = ptr_array.get();

    EdgeLoader loader { handle };
    uint64_t num_edges; uint64_t i = 0;
    uint64_t num_blocks = 8;
    while( (num_edges = loader.load(array, array_sz) ) > 0 ){
        if(i +1 < num_blocks){
            ASSERT_EQ(num_edges, 16ull);
        } else {
            ASSERT_EQ(num_edges, 8ull);
        }
        i++;
    }

    ASSERT_EQ(i, num_blocks);

    handle.close();
}

static void validate_edge(EdgeReader& reader, uint64_t expected_source_id, uint64_t expected_destination_id, double expected_weight){
    graph::WeightedEdge edge;
    bool has_read_edge = reader.read_edge(edge);
    ASSERT_TRUE(has_read_edge);
    ASSERT_EQ(edge.source(), expected_source_id);
    ASSERT_EQ(edge.destination(), expected_destination_id);
    ASSERT_EQ(edge.weight(), expected_weight);
}

TEST(Graphlog, EdgeReader){
    EdgeReader reader { path_graph };

    /**
     * Adapt the script below to generate the following list:
     *
     * #!/usr/bin/env perl
     *
     * use strict;
     * use warnings;
     *
     * open(my $f, "input.txt"); # create a log (input.txt) in the generator of the form [src: <id>, dst: <id>, weight: <value>]
     * my $matches = 0;
     *
     * while(my $line = <$f>){
     *     if($line =~ /\[src: (\d+), dst: (\d+), weight: (-?[\d\.]+)\]/){
     *         print("validate_edge(reader, $1, $2, $3);\n");
     *         $matches ++;
     *     }
     * }
     *
     * close($f);
     *
     * print("matches: $matches\n");
     *
     */

    validate_edge(reader, 2, 3, 0);
    validate_edge(reader, 5, 8, 0.12);
    validate_edge(reader, 6, 8, 0);
    validate_edge(reader, 4, 8, 0);
    validate_edge(reader, 2, 7, 0);
    validate_edge(reader, 6, 10, 0);
    validate_edge(reader, 3, 6, 0);
    validate_edge(reader, 3, 5, 0);
    validate_edge(reader, 3, 11, 0);
    validate_edge(reader, 4, 6, 0);
    validate_edge(reader, 5, 6, 0);
    validate_edge(reader, 7, 9, 0.36);
    validate_edge(reader, 4, 6, -1);
    validate_edge(reader, 9, 11, 0);
    validate_edge(reader, 5, 6, -1);
    validate_edge(reader, 10, 11, 0);
    validate_edge(reader, 6, 10, -1);
    validate_edge(reader, 1, 3, 0);
    validate_edge(reader, 3, 6, -1);
    validate_edge(reader, 1, 6, 0);
    validate_edge(reader, 3, 5, -1);
    validate_edge(reader, 3, 5, 0.5);
    validate_edge(reader, 6, 8, -1);
    validate_edge(reader, 5, 6, 0);
    validate_edge(reader, 2, 3, -1);
    validate_edge(reader, 4, 6, 0);
    validate_edge(reader, 1, 6, -1);
    validate_edge(reader, 8, 10, 0);
    validate_edge(reader, 8, 10, -1);
    validate_edge(reader, 6, 8, 0);
    validate_edge(reader, 10, 11, -1);
    validate_edge(reader, 3, 8, 0.32);
    validate_edge(reader, 9, 11, -1);
    validate_edge(reader, 2, 8, 0);
    validate_edge(reader, 5, 6, -1);
    validate_edge(reader, 3, 6, 0);
    validate_edge(reader, 6, 8, -1);
    validate_edge(reader, 6, 8, 0);
    validate_edge(reader, 3, 11, -1);
    validate_edge(reader, 9, 10, 0);
    validate_edge(reader, 2, 7, -1);
    validate_edge(reader, 6, 7, 0.53);
    validate_edge(reader, 3, 6, -1);
    validate_edge(reader, 3, 6, 0);
    validate_edge(reader, 2, 8, -1);
    validate_edge(reader, 3, 4, 0);
    validate_edge(reader, 9, 10, -1);
    validate_edge(reader, 2, 4, 0);
    validate_edge(reader, 4, 8, -1);
    validate_edge(reader, 8, 11, 0);
    validate_edge(reader, 6, 8, -1);
    validate_edge(reader, 6, 8, 0.64);
    validate_edge(reader, 4, 6, -1);
    validate_edge(reader, 2, 6, 0);
    validate_edge(reader, 1, 3, -1);
    validate_edge(reader, 2, 5, 0);
    validate_edge(reader, 2, 4, -1);
    validate_edge(reader, 4, 7, 0);
    validate_edge(reader, 2, 6, -1);
    validate_edge(reader, 1, 7, 0);
    validate_edge(reader, 4, 7, -1);
    validate_edge(reader, 2, 3, 0.9);
    validate_edge(reader, 1, 7, -1);
    validate_edge(reader, 5, 9, 0);
    validate_edge(reader, 3, 4, -1);
    validate_edge(reader, 7, 8, 0);
    validate_edge(reader, 5, 9, -1);
    validate_edge(reader, 1, 9, 0);
    validate_edge(reader, 7, 8, -1);
    validate_edge(reader, 1, 11, 0);
    validate_edge(reader, 3, 6, -1);
    validate_edge(reader, 2, 4, 0.69);
    validate_edge(reader, 1, 11, -1);
    validate_edge(reader, 6, 9, 0);
    validate_edge(reader, 2, 5, -1);
    validate_edge(reader, 5, 6, 0);
    validate_edge(reader, 6, 9, -1);
    validate_edge(reader, 3, 9, 0);
    validate_edge(reader, 1, 9, -1);
    validate_edge(reader, 2, 9, 0);
    validate_edge(reader, 3, 9, -1);
    validate_edge(reader, 5, 6, -1);
    validate_edge(reader, 5, 6, 0.63);
    validate_edge(reader, 4, 11, 0);
    validate_edge(reader, 2, 9, -1);
    validate_edge(reader, 9, 11, 0);
    validate_edge(reader, 4, 11, -1);
    validate_edge(reader, 4, 7, 0);
    validate_edge(reader, 8, 11, -1);
    validate_edge(reader, 4, 8, 0);
    validate_edge(reader, 4, 7, -1);
    validate_edge(reader, 3, 4, 0.13);
    validate_edge(reader, 4, 8, -1);
    validate_edge(reader, 2, 10, 0);
    validate_edge(reader, 9, 11, -1);
    validate_edge(reader, 3, 10, 0);
    validate_edge(reader, 3, 10, -1);
    validate_edge(reader, 6, 9, 0);
    validate_edge(reader, 6, 9, -1);
    validate_edge(reader, 1, 6, 0);
    validate_edge(reader, 1, 6, -1);
    validate_edge(reader, 6, 10, 0.63);
    validate_edge(reader, 2, 10, -1);
    validate_edge(reader, 3, 11, 0);
    validate_edge(reader, 3, 11, -1);
    validate_edge(reader, 5, 7, 0);
    validate_edge(reader, 5, 7, -1);
    validate_edge(reader, 6, 9, 0);
    validate_edge(reader, 6, 9, -1);
    validate_edge(reader, 7, 8, 0);
    validate_edge(reader, 7, 8, -1);
    validate_edge(reader, 6, 9, 0.23);
    validate_edge(reader, 3, 10, 0);
    validate_edge(reader, 3, 10, -1);
    validate_edge(reader, 2, 10, 0);
    validate_edge(reader, 2, 10, -1);
    validate_edge(reader, 1, 6, 0);
    validate_edge(reader, 1, 6, -1);
    validate_edge(reader, 3, 9, 0);
    validate_edge(reader, 3, 9, -1);

    { // stream depleted
        graph::WeightedEdge edge;
        ASSERT_FALSE( reader.read_edge(edge) );
        ASSERT_FALSE( reader.read_edge(edge) );
    }
}
