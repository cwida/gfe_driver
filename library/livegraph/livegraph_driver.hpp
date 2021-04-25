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

#include <atomic>
#include <chrono>
#include "library/interface.hpp"

namespace gfe::library {

/**
 * Wrapper to evaluate the LiveGraph library
 */
class LiveGraphDriver : public virtual UpdateInterface, public virtual GraphalyticsInterface {
    LiveGraphDriver(const LiveGraphDriver&) = delete;
    LiveGraphDriver& operator=(const LiveGraphDriver&) = delete;

protected:
    void* m_pImpl; // pointer to the LiveGraph handle
    void* m_pHashMap; // pointer to the TBB HashMap to translate the vertex identifiers into the dense IDs for livegraph
    const bool m_is_directed; // whether the underlying graph is directed or undirected
    const bool m_read_only; // whether to used read only transactions for graphalytics
    std::atomic<uint64_t> m_num_vertices {0}; // keep track of the total number of vertices
    std::atomic<uint64_t> m_num_edges {0}; // keep track of the total number fo edges
    std::chrono::seconds m_timeout {0}; // the budget to complete each of the algorithms in the Graphalytics suite


    // Retrieve the internal vertex ID for the given external vertex. If the vertex does not exist, it raises an internal error
    uint64_t ext2int(uint64_t external_vertex_id) const;

    // Retrieve the internal vertex ID for the given internal vertex ID. If the vertex does not exist, it returns uint64_t::max()
    uint64_t int2ext(void* transaction, uint64_t internal_vertex_id) const;

    // Helper for Graphalytics: translate the logical IDs into external IDs
    template <typename T>
    std::vector<std::pair<uint64_t, T>> translate(void* /* transaction object */ lgtxn, const T* __restrict data, uint64_t data_sz);

    // Helper, save the content of the vector to the given output file
    template <typename T, bool negative_scores = true>
    void save_results(const std::vector<std::pair<uint64_t, T>>& result, const char* dump2file);
public:
    /**
     * Create an instance of LiveGraph
     * @param is_directed: whether the underlying graph should be directed or undirected
     * @param read_only: whether to use read-only transactions for the algorithms in Graphalytics
     */
    LiveGraphDriver(bool is_directed, bool read_only = true);

    /**
     * Destructor
     */
    virtual ~LiveGraphDriver();

    /**
     * Get the number of edges contained in the graph
     */
    virtual uint64_t num_edges() const;

    /**
     * Get the number of nodes stored in the graph
     */
    virtual uint64_t num_vertices() const;

    /**
     * Returns true if the given vertex is present, false otherwise
     */
    virtual bool has_vertex(uint64_t vertex_id) const;

    /**
     * Returns the weight of the given edge is the edge is present, or NaN otherwise
     */
    virtual double get_weight(uint64_t source, uint64_t destination) const;

    /**
     * Check whether the graph is directed
     */
    virtual bool is_directed() const;

    /**
     * Impose a timeout on each graph computation. A computation that does not terminate by the given seconds will raise a TimeoutError.
     */
    virtual void set_timeout(uint64_t seconds);

    /**
     * Add the given vertex to the graph
     * @return true if the vertex has been inserted, false otherwise (that is, the vertex already exists)
     */
    virtual bool add_vertex(uint64_t vertex_id);

    /**
     * Remove the mapping for a given vertex. The actual internal vertex is not removed from the adjacency list.
     * @param vertex_id the vertex to remove
     * @return true if a mapping for that vertex existed, false otherwise
     */
    virtual bool remove_vertex(uint64_t vertex_id);

    /**
     * Add the given edge in the graph if it doesn't exist
     * @return true if the edge has been inserted, false if this edge already exists or one of the referred
     *         vertices does not exist.
     */
    virtual bool add_edge(gfe::graph::WeightedEdge e);

    /**
     * Add the given edge in the graph. Implicitly create the referred vertices if they do not already exist.
     * If the edge already exists, its weight is updated.
     * @return always true.
     */
    virtual bool add_edge_v2(gfe::graph::WeightedEdge e);

    /**
     * Remove the given edge from the graph
     * @return true if the given edge has been removed, false otherwise (e.g. this edge does not exist)
     */
    virtual bool remove_edge(gfe::graph::Edge e);

    /**
     * Dump the content of the graph to given stream.
     */
    virtual void dump_ostream(std::ostream& out) const;

    /**
     * Retrieve the opaque objects for the internal LiveGraph & Vertex Dictionary handles.
     * For Debugging & Testing only
     */
    void* livegraph(); // lg::Graph*
    void* vertex_dictionary(); // tbb::concurrent_hash_map<uint64_t, lg::vertex_t>

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

} // namespace
