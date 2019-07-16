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

#include "common/error.hpp"
#include "library/interface.hpp"

#include <cinttypes>
#include <shared_mutex>
#include <unordered_map>
#include <vector>

namespace library {

// Generic exception thrown by this class
DEFINE_EXCEPTION(AdjacencyListError);

/**
 * Sequential and base implementation of the interface, for testing purposes.
 * The class is thread-safe, but all exposed operations are serialised and sequential.
 */
class AdjacencyList : public virtual UpdateInterface, public virtual LoaderInterface, public virtual GraphalyticsInterface {
    using EdgeList = std::vector</* edge */ std::pair< /* destination */ uint64_t,  /* weight */ double>>;
    using EdgePair = std::pair< /* outgoing edges */ EdgeList, /* incoming edges */ EdgeList >;
    using NodeList = std::unordered_map</* src vertex */ uint64_t, /* outgoing & incoming edges */ EdgePair>;
    NodeList m_adjacency_list;
    uint64_t m_num_edges = 0; // number of directed edges
    const bool m_is_directed; // whether the graph is directed or not

    using mutex_t = std::shared_mutex;
    mutable mutex_t m_mutex; // read-write mutex

    // Get the list of incoming edges
    const EdgeList& get_incoming_edges(uint64_t vertex_id) const;
    const EdgeList& get_incoming_edges(const NodeList::iterator& it) const;
    const EdgeList& get_incoming_edges(const NodeList::const_iterator& it) const;
    const EdgeList& get_incoming_edges(const NodeList::mapped_type& adjlist_value) const;

    // Get the degree of the given vertex v, that is | in(v) U out(v) | with U = set union
    uint64_t get_degree(uint64_t vertex_id) const;

    // Apply the function F(vertex) to all outgoing & incoming edges
    template<typename Function>
    void for_all_edges(const EdgePair&, Function F);

    // Check if a directed edge between vertex1 and vertex2 exists
    bool check_directed_edge_exists(uint64_t vertex1, uint64_t vertex2);

    // Assume the lock has already been acquired
    bool add_vertex0(uint64_t vertex_id);
    bool delete_vertex0(uint64_t vertex_id);
    bool add_edge0(graph::WeightedEdge e);
    bool delete_edge0(graph::Edge e);


public:
    /**
     * Initialise the graph instance
     */
    AdjacencyList(bool is_directed);

    /**
     * Destructor
     */
    ~AdjacencyList();

    /**
     * Get the number of edges contained in the graph
     */
    virtual uint64_t num_edges() const;

    /**
     * Get the number of nodes stored in the graph
     */
    virtual uint64_t num_vertices() const;

    /**
     * Is the graph directed?
     */
    virtual bool is_directed() const;

    /**
     * Returns true if the given vertex is present, false otherwise
     */
    virtual bool has_vertex(uint64_t vertex_id) const;

    /**
     * Retrieve the weight associated to the given edge, or -1 if the given edge does not exist
     */
    virtual double get_weight(uint64_t source, uint64_t destination) const;

    /**
     * Dump the content of the graph to the given output stream
     */
    virtual void dump_ostream(std::ostream& out) const;

    /**
     * Add the given vertex to the graph
     * @return true if the vertex has been inserted, false otherwise
     */
    virtual bool add_vertex(uint64_t vertex_id);

    /**
     * Remove the given vertex and all edges attached to it.
     * @return true in case of success, false otherwise
     */
    virtual bool delete_vertex(uint64_t vertex_id);

    /**
     * Add the given edge in the graph
     * @return true if the edge has been inserted or updated, false in case of error
     */
    virtual bool add_edge(graph::WeightedEdge e);

    /**
     * Remove the given edge from the graph
     * @return true if the given edge has been removed, false otherwise (e.g. this edge does not exist)
     */
    virtual bool delete_edge(graph::Edge e);

    /**
     * Load the whole graph representation from the given path
     */
    virtual void load(const std::string& path);

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

