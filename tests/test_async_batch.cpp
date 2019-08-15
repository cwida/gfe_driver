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

#include <cstdlib> // getenv
#include <memory>
#include <unordered_set>

#include "experiment/details/async_batch.hpp"
#include "graph/edge.hpp"
#include "library/baseline/adjacency_list.hpp"

using namespace experiment::details;
using namespace graph;
using namespace library;
using namespace std;

TEST(AsyncBatch, Sanity){
    auto library = make_shared<AdjacencyList>(/* directed = */ true);
    AsyncBatch batch { library.get(), 1, /* num batches */ 3, /* batch_sz */ 4};

    for(uint64_t i = 1; i <= 8; i++){
        library->add_vertex(i);
    }

    for(uint64_t i = 1; i < 8; i++){
        for(uint64_t j = i +1; j <= 8; j++){
            batch.add_edge(graph::WeightedEdge{i, j, (double) i + j});
        }
    }
    batch.flush(true);

//    library->dump();

    for(uint64_t i = 1; i < 8; i++){
        for(uint64_t j = i +1; j <= 8; j++){
            ASSERT_TRUE(library->has_edge(i, j));
            ASSERT_EQ(library->get_weight(i, j), i +j);
        }
    }
}
