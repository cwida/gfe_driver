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
#include <cstddef>
#include <cstdint>

#include "common/spinlock.hpp"
#include "library/interface.hpp"

namespace gfe::library {

/**
 * Wrapper to evaluate the GraphOne library
 */
class GraphOne : public virtual UpdateInterface, public virtual GraphalyticsInterface {
    GraphOne(const GraphOne&) = delete;
    GraphOne& operator=(const GraphOne&) = delete;

    // you might wonder where it the handle to the actual implementation? Well, it's a global
    // variable, as mandated by the library, named g and defined in base.cpp.
    // ...

    const bool m_is_directed; // whether the underlying graph is directed or undirected
    const bool m_translate_vertex_ids; // whether to use the mapping from external to internal for the vertex identifiers
    const bool m_blind_writes; // whether we check the presence of prior existing elements before performing an insertion
    std::chrono::seconds m_timeout { 0 }; // the budget to complete each of the algorithms in the Graphalytics suite
    std::atomic<uint64_t> m_num_vertices { 0 }; // total number of vertices in the graph
    uint64_t m_num_edges { 0 }; // total number of edges in the graph (not just those archived)
    struct PaddedLock { // to avoid false sharing
        common::SpinLock m_lock;
        uint64_t padding[7];
    };
    common::SpinLock m_mutex_vtx; // mutex for the vertex dictionary
    uint64_t m_num_edge_locks { 0 }; // the size of the array m_edge_locks;
    PaddedLock* m_edge_locks { nullptr }; // battery of lock, used for the edge updates, when blink_writes is disabled

    // Protect the single invocation to #do_update
    common::SpinLock m_do_update_lock; // The interface from GraphOne for edges updates is not thread safe

    // This is an optimisation. The implementation of #get_weight() is rather slow if it has to create a static view every time
    // it is invoked. We keep a cache of the last static view created
//    void* m_ptr_static_view_cached { nullptr }; // check whether the last created static view
//    std::atomic<bool> m_static_view_is_valid = false; // whether the last created static view is still valid

    // Insert/remove an edge into GraphOne. The weight is ignored if the operation is a deletion
    void do_update(bool is_insert, uint64_t logical_source_id, uint64_t logical_destination_id, double weight = 0.0);

    /**
     * Find the given edge in both the write and read store. This method does not acquire the m_edge_locks, this is supposed
     * to be done by the calling function
     * @param src the source vertex of the edge, as internal logical id (sid_t)
     * @param dst the destination vertex of the edge, as internal logical id (sid_t)
     * @param out_found output value, whether a record with the given source and destination has been found, regardless if it's an actual insertion or a tombstone for a deletion
     * @param out_weight output value, if the edge exists, the weight associated to the given edge
     * @return true if the edge exists in the storage, false otherwise
     */
    bool find_edge(uint64_t src, uint64_t dst, bool* out_record_found = nullptr, double* out_weight = nullptr) const;

    // Internal implementation of the method #get_weight
    double get_weight_impl(uint64_t v0, uint64_t v1) const;

    // Given an external vertex ID, retrieve the internal (logical) vertex ID, or raise an exception if the mapping does not exist
    uint64_t vtx_ext2int(uint64_t external_vertex_id) const;

public:

    /**
     * Create an instance of the wrapper for GraphOne
     * @param is_graph_directed: whether the underlying graph should be directed or undirected
     * @param use_vertex2id_mapping: true for sparse graphs, false for dense graphs. A dense graph assumes the vertex IDs are in
     *        a contiguous sequence in [0, |V|), with |V| the total number of vertices inserted so far, while in a sparse graph
     *        a vertex ID can be any positive integer.
     * @param blind_writes: it affects consistency of edge insertions/deletions. If false, every time an insertion is performed,
     *        it checks whether that edge already exists, and, similarly, it only allows edge deletions for edges that already
     *        exist. This kind of explicit consistency is rather expensive on the current version of GraphOne. Enabling blind
     *        writes assumes that all edge insertions and deletions are correct, i.e. always refer to non existing or existing edges
     *        respectively, and the methods #add_edge and #remove_add always report success, unless the related vertices do not
     *        exist.
     * @param max_num_vertices: in GraphOne the array of vertices must be statically allocated in advance, with a fixed size.
     *        In a server environment, we set the capacity to 4G (= 2^32). For a workstation, this value is too high and does
     *        not enable creating an instance for testing purposes. This parameter allows to explicitly set the capacity of
     *        the vertex array.
     */
    GraphOne(bool is_graph_directed, bool use_vertex2id_mapping, bool blind_writes, uint64_t max_num_vertices);

    /**
     * Destructor (dummy)
     */
    virtual ~GraphOne();

    /**
     * Dump the content of the graph to given stream
     */
    void dump_ostream(std::ostream& out) const;

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
     * Whether blind writes are enabled or not. Blind writes refer to whether each edge insertion or deletion check
     * its prior presence in the storage, to ensure consistency. If blind writes are enabled, consistency cannot be
     * enforced, and, for instance, removing an edge that does not already exist may cause a crash.
     */
    bool has_blind_writes() const;

    /**
     * Returns the weight of the given edge is the edge is present, or NaN otherwise
     */
    virtual double get_weight(uint64_t source, uint64_t destination) const;

    /**
     * Check whether the graph is directed
     */
    virtual bool is_directed() const;

    /**
     * Check whether the graph is undirected
     */
    virtual bool is_undirected() const;

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
     * Add the given edge in the graph. The implementation does not check whether this edge already exists,
     * adding a new edge always.
     * @return always true when both the source & the destination vertices already exist, false otherwise
     */
    virtual bool add_edge(gfe::graph::WeightedEdge e);

    /**
     * Remove the given edge from the graph. There is no way to check whether the operation actually succeeded
     * in this implementation of GraphOne. Attempting to remove an edge that does not exist may result in a crash.
     * @return always true when both the source & the destination vertices already exist, false otherwise
     */
    virtual bool remove_edge(gfe::graph::Edge e);

    /**
     * Flush all vertices from the buffer into the read-only adjacency lists
     */
    virtual void build();

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
