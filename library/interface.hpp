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

#include <iostream>
#include <memory>
#include <ostream>
#include <vector>

#include "graph/edge.hpp"

namespace library {

// Forward declarations
class Interface;
class UpdateInterface;
class LoaderInterface;

/**
 * The list of registered implementations
 */
class ImplementationManifest {
public:
    std::string m_name; // unique id
    std::string m_description; // description of the instance, showed in the help screen ( -h )
    std::unique_ptr<Interface> (*m_factory)(void); // factory method to generate an instance of this implementation

    ImplementationManifest(const std::string& name, const std::string& description, std::unique_ptr<Interface> (*factory)(void));
};
std::vector<ImplementationManifest> implementations();

/**
 * Base interface, implemented by all systems
 */
class Interface {
public:
    /**
     * Dummy constructor
     */
    Interface();

    /**
     * Virtual destructor
     */
    virtual ~Interface();

    /**
     * Dump the content of the graph to stdout
     */
    virtual void dump() const = 0;

    /**
     * Thread initialisation callbacks
     */
    // Invoked at the start of the experiment by the controller thread with the number of threads that will be used
    virtual void on_main_init(int num_threads);
    // Invoked by each worker thread separately, with a unique thread id, in [0, num_threads)
    virtual void on_thread_init(int thread_id);
    // Invoked by each worker thread separately, with the same thread id given at on_thread_init
    virtual void on_thread_destroy(int thread_id);
    // Invoked at the end of the experiment by the controller thread
    virtual void on_main_destroy();

    /**
     * Get the number of edges contained in the graph
     */
    virtual uint64_t num_edges() const = 0;

    /**
     * Get the number of nodes stored in the graph
     */
    virtual uint64_t num_vertices() const = 0;

    /**
     * Returns true if the given vertex is present, false otherwise
     */
    virtual bool has_vertex(uint64_t vertex_id) const = 0;

    /**
     * Returns true if the given edge is present, false otherwise
     */
    virtual bool has_edge(uint64_t source, uint64_t destination) const;

    /**
     * Returns the weight of the given edge is the edge is present, or -1 otherwise
     */
    virtual int64_t get_weight(uint64_t source, uint64_t destination) const = 0;
};


/**
 * Update interface
 */
class UpdateInterface : public virtual Interface {
public:
    /**
     * Add the given vertex to the graph
     * @return true if the edge has been inserted, false otherwise (e.g. the vertex already exists)
     */
    virtual bool add_vertex(uint64_t vertex_id) = 0;

    /**
     * Remove the given vertex from the graph. The implementation may be unable to remove the vertex
     * if there are still edges attached or it may go ahead the vertex and all edges attached. This
     * behaviour is implementation dependent.
     * @param vertex_id the vertex to remove
     * @return true in case of success, false otherwise.
     */
    virtual bool delete_vertex(uint64_t vertex_id) = 0;

    /**
     * Add the given edge in the graph
     * @return true if the edge has been inserted, false otherwise (e.g. this edge already exists)
     */
    virtual bool add_edge(graph::WeightedEdge e) = 0;

    /**
     * Remove the given edge from the graph
     * @return true if the given edge has been removed, false otherwise (e.g. this edge does not exist)
     */
    virtual bool delete_edge(graph::Edge e) = 0;
};


/**
 * Load the graph from a file in the disk
 */
class LoaderInterface : public virtual Interface {
public:
    /**
     * Load the whole graph representation from the given path
     */
    virtual void load(const std::string& path) = 0;
};


/**
 * The six algorithms required by the Graphalytics benchmark suite
 * See https://github.com/ldbc/ldbc_graphalytics_docs/
 */
class GraphalyticsInterface : public virtual Interface {
public:
    /**
     * Perform a BFS from source_vertex_id to all the other vertices in the graph.
     * @param source_vertex_id the vertex where to start the search
     * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
     */
    virtual void bfs(uint64_t source_vertex_id, const char* dump2file = nullptr) = 0;

    /**
     * Execute the PageRank algorithm for the specified number of iterations.
     *
     * @param num_iterations the number of iterations to execute the algorithm
     * @param damping_factor weight for the PageRank algorithm, it affects the score associated to the sink nodes in the graphs
     * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
     */
    virtual void pagerank(uint64_t num_iterations, double damping_factor = 0.85, const char* dump2file = nullptr) = 0;

    /**
     * Weakly connected components (WCC), associate each node to a connected component of the graph
     * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
     */
    virtual void wcc(const char* dump2file = nullptr) = 0;

    /**
     * Community Detection using Label-Propagation. Associate a label to each vertex of the graph, according to its neighbours.
     * @param max_iterations maximum number of iterations to perform
     * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
     */
    virtual void cdlp(uint64_t max_iterations, const char* dump2file = nullptr) = 0;

    /**
     * Local clustering coefficient. Associate to each vertex the ratio between the number of its outgoing edges and the number of
     * possible remaining edges.
     * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
     */
    virtual void lcc(const char* dump2file = nullptr) = 0;

    /**
     * Single-source shortest paths. Compute the weight related to the shortest path from the source to any other vertex in the graph.
     * @param source_vertex_id the vertex where to start the search
     * @param dump2file if not null, dump the result in the given path, following the format expected by the benchmark specification
     */
    virtual void sssp(uint64_t source_vertex_id, const char* dump2file = nullptr) = 0;
};

} // namespace library

