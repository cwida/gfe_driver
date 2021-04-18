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

#pragma once

#include "teseo_driver.hpp"

namespace gfe::library {

/**
 * Implementation of the Graphalytics kernels assuming that the vertices in the graphs have been
 * relabelled in advance in [0, N) (dense graphs). The algorithm implementations are the same
 * of the TeseoDriver, except that the iterator scans the real vertices, rather than the logical IDs.
 */
class TeseoRealVertices : public TeseoDriver {
protected:
    // Helper, save the content of the vector to the given output file
    template <typename T>
    void save_results(const T* __restrict result, uint64_t result_sz, const char* dump2file);

public:
    /**
     * Create an instance of Teseo
     * @param is_directed: whether the underlying graph should be directed or undirected
     * @param read_only: whether to use read-only transactions for the algorithms in graphalytics
     */
    TeseoRealVertices(bool is_directed, bool read_only = true);

    /**
     * Perform a BFS from source_vertex_id to all the other vertices in the graph.
     * @param source_vertex_id the vertex where to start the search
     * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
     */
    virtual void bfs(uint64_t source_vertex_id, const char* dump2file = nullptr);

    /**
     * Execute the PageRank algorithm for the specified number of iterations.
     *
     * @param num_iterations the number of iterations to execute the algorithm
     * @param damping_factor weight for the PageRank algorithm, it affects the score associated to the sink nodes in the graphs
     * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
     */
    virtual void pagerank(uint64_t num_iterations, double damping_factor = 0.85, const char* dump2file = nullptr);

    /**
     * Weakly connected components (WCC), associate each node to a connected component of the graph
     * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
     */
    virtual void wcc(const char* dump2file = nullptr);

    /**
     * Community Detection using Label-Propagation. Associate a label to each vertex of the graph, according to its neighbours.
     * @param max_iterations max number of iterations to perform
     * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
     */
    virtual void cdlp(uint64_t max_iterations, const char* dump2file = nullptr);

    /**
     * Local clustering coefficient. Associate to each vertex the ratio between the number of its outgoing edges and the number of
     * possible remaining edges.
     * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
     */
    virtual void lcc(const char* dump2file = nullptr);

    /**
     * Single-source shortest paths. Compute the weight related to the shortest path from the source to any other vertex in the graph.
     * @param source_vertex_id the vertex where to start the search
     * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
     */
    virtual void sssp(uint64_t source_vertex_id, const char* dump2file = nullptr);
};

/**
 * Specialised implementation of the LCC kernel
 */
class TeseoRealVerticesLCC : public TeseoRealVertices {
    TeseoRealVerticesLCC(const TeseoRealVerticesLCC& ) = delete;
    TeseoRealVerticesLCC& operator=(const TeseoRealVerticesLCC& ) = delete;

public:
    /**
     * Constructor
     * @param is_directed whether the graph is directed
     * @param read_only whether to create the transactions for the algorithms in Graphalytics as read-only
     */
    TeseoRealVerticesLCC(bool is_directed, bool read_only = true);

    /**
     * Specialised implementation of the kernel LCC
     */
    virtual void lcc(const char* dump2file = nullptr);
};

} // namespace
