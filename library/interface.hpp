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

#include "common/error.hpp"
#include "graph/edge.hpp"

namespace gfe::library {

// Forward declarations
class Interface;
class UpdateInterface;
class LoaderInterface;

// Raised if the computation did not complete in the given time budget.
DEFINE_EXCEPTION(TimeoutError);

/**
 * The manifest associated the name of a system implementation, that can be evaluated, with a factory method to create an instance
 * of the given system. A manifest is an entry in the list of all available implementations, that can be retrieved using the function
 * gfe::library::implementations().
 */
class ImplementationManifest {
public:
    std::string m_name; // unique id
    std::string m_description; // description of the instance, showed in the help screen ( -h )
    std::unique_ptr<Interface> (*m_factory)(bool is_graph_directed); // factory method to generate an instance of this implementation

    /**
     * Creates a new manifest, that is an association between the name of a system and a factory function to create an instance of
     * the related system.
     * @param name: the name of the system, as invoked externally by command line and saved in the final database with the results
     * @param description: a short description of the system, showed in the command line when the program is invoked with --help
     * @param factory: the factory function, to actually an instance of the system, either for directed or undirected graphs
     */
    ImplementationManifest(const std::string& name, const std::string& description, std::unique_ptr<Interface> (*factory)(bool is_graph_directed));
};

// Retrieve the list of all implementations that can be evaluated, together with the factory methods to create them
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
     * Dump the content of the graph to given stream
     */
    virtual void dump_ostream(std::ostream& out) const = 0;

    /**
     * Dump the content of the graph to stdout
     */
    virtual void dump() const;

    /**
     * Dump the content of the graph to the given path
     */
    virtual void dump(const std::string& path) const;
    virtual void dump(const char* path) const; // helps jitting from the debugger

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
     * Returns the weight of the given edge is the edge is present, or NaN otherwise
     */
    virtual double get_weight(uint64_t source, uint64_t destination) const = 0;

    /**
     * Check whether the graph is directed
     */
    virtual bool is_directed() const = 0;

    /**
     * Check whether the graph is undirected
     */
    virtual bool is_undirected() const;

    /**
     * Impose a timeout on each graph computation. A computation that does not terminate by the given seconds will raise a TimeoutError.
     */
    virtual void set_timeout(uint64_t seconds) = 0;

    // To assess the overhead of the compaction phase in LLAMA
    virtual void updates_start();
    virtual void updates_stop();

    /**
     * Check whether we are allowed to validate the updates performed
     */
    virtual bool can_be_validated() const;
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
 * Retrieve a random vertex ID
 */
class RandomVertexInterface : public virtual Interface {
public:
    /**
     * Get a random vertex ID
     */
    virtual uint64_t get_random_vertex_id() const = 0;
};

/**
 * Update interface
 */
class UpdateInterface : public virtual Interface, public virtual LoaderInterface {
private:
    /**
     * Helper function for the implementation of bool batch(...);
     * Constantly invoke action(edge) until either it returns true or a timeout has expired. If timeout
     * expired, the function throws a TimeoutError.
     */
    template<typename Action, typename Edge>
    void batch_try_again(Action action, Edge edge);

public:
    /**
     * Add the given vertex to the graph
     * @return true if the vertex has been inserted, false otherwise (that is, the vertex already exists)
     */
    virtual bool add_vertex(uint64_t vertex_id) = 0;

    /**
     * Remove the given vertex from the graph. The implementation may be unable to remove the vertex
     * if there are still edges attached or it may go ahead the vertex and all edges attached. This
     * behaviour is implementation dependent.
     * @param vertex_id the vertex to remove
     * @return true in case of success, false otherwise.
     */
    virtual bool remove_vertex(uint64_t vertex_id) = 0;

    /**
     * Add the given edge in the graph
     * @return true if the edge has been inserted, false if this edge already exists or one of the referred
     *    vertices does not exist.
     */
    virtual bool add_edge(gfe::graph::WeightedEdge e) = 0;

    /**
     * Add the given edge in the graph. Implicitly create the referred vertices if they do not already exist
     * @return true if the edge has been inserted, false otherwise (e.g. this edge already exists)
     */
    virtual bool add_edge_v2(gfe::graph::WeightedEdge e) = 0;

    /**
     * Remove the given edge from the graph
     * @return true if the given edge has been removed, false otherwise (e.g. this edge does not exist)
     */
    virtual bool remove_edge(gfe::graph::Edge e) = 0;

    /**
     * Load the whole graph representation from the given path
     */
    virtual void load(const std::string& path) override; // default implementation provided in terms of #add_vertex and #add_edge

    /**
     * Create a new snapshot. By default this operation is a `nop'.
     * In LLAMA, it creates a new level, moving all pending updates in the write store into a new delta in the read-only store.
     */
    virtual void build(); // default impl. => `nop'

    /**
     * Returns the total number of snaphots/levels/deltas created in the LSM/delta based implementation.
     */
    virtual uint64_t num_levels() const;

    /**
     * Perform a batch of edge insertions/deletions.
     * -- LIBRARY IMPLEMENTATIONS SHALL NOT OVERRIDE THIS METHOD: this is only used by the driver in client-server
     * mode. The point is not to measure the performance of ``batch updates'' in the library, but to amortize the
     * cost of many RPC calls over the network.
     *
     * @param array the list of edge updates, insertions/deletions
     * @param array_sz the size of the list of edge updates
     * @param force if true, expect to perform all updates successfully and wait indefinitely until each operation is performed (due to async threads)
     * @return true if all updates have been performed, and false if one of them failed
     */
    struct SingleUpdate {
        uint64_t m_source; // the source vertex
        uint64_t m_destination; // the destination vertex
        double m_weight; // if < 0, this is an edge removal, otherwise it's an edge insertion with the given weight
    };
    virtual bool batch(const SingleUpdate* array, size_t array_sz, bool force = true);
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
     * @param max_iterations max number of iterations to perform
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

} // namespace

