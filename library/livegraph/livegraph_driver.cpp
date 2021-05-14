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

#include "livegraph_driver.hpp"

#include <atomic>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <mutex>
#include <limits>
#include <sstream>
#include <thread>
#include <unordered_set>

#include "common/system.hpp"
#include "common/timer.hpp"
#include "tbb/concurrent_hash_map.h"
#include "third-party/gapbs/gapbs.hpp"
#include "third-party/libcuckoo/cuckoohash_map.hh"
#include "third-party/livegraph/livegraph.hpp"
#include "utility/timeout_service.hpp"

using namespace common;
using namespace libcuckoo;
using namespace std;

#define LiveGraph reinterpret_cast<lg::Graph*>(m_pImpl)
using vertex_dictionary_t = tbb::concurrent_hash_map<uint64_t, lg::vertex_t>;
#define VertexDictionary reinterpret_cast<vertex_dictionary_t*>(m_pHashMap)

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG
namespace gfe { extern mutex _log_mutex [[maybe_unused]]; }
#define COUT_DEBUG_FORCE(msg) { std::scoped_lock<std::mutex> lock{::gfe::_log_mutex}; std::cout << "[LiveGraphDriver::" << __FUNCTION__ << "] [Thread #" << common::concurrency::get_thread_id() << "] " << msg << std::endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif

namespace gfe::library {

/*****************************************************************************
 *                                                                           *
 *  Init                                                                     *
 *                                                                           *
 *****************************************************************************/
LiveGraphDriver::LiveGraphDriver(bool is_directed, bool read_only) : m_pImpl(nullptr), m_pHashMap(nullptr), m_is_directed(is_directed), m_read_only(read_only) {
    m_pImpl = new lg::Graph();
    m_pHashMap = new tbb::concurrent_hash_map<uint64_t, /* vertex_t */ uint64_t>();
}

LiveGraphDriver::~LiveGraphDriver(){
    delete LiveGraph; m_pImpl = nullptr;
    delete VertexDictionary; m_pHashMap = nullptr;
}

/*****************************************************************************
 *                                                                           *
 *  Properties                                                               *
 *                                                                           *
 *****************************************************************************/
bool LiveGraphDriver::is_directed() const {
    return m_is_directed;
}

uint64_t LiveGraphDriver::num_edges() const {
    return m_num_edges;
}

uint64_t LiveGraphDriver::num_vertices() const {
    return m_num_vertices;
}

void LiveGraphDriver::set_timeout(uint64_t seconds) {
    m_timeout = chrono::seconds{ seconds };
}

void* LiveGraphDriver::livegraph() {
    return m_pImpl;
}

void* LiveGraphDriver::vertex_dictionary() {
    return m_pHashMap;
}

uint64_t LiveGraphDriver::ext2int(uint64_t external_vertex_id) const {
    vertex_dictionary_t::const_accessor accessor;
    if ( VertexDictionary->find(accessor, external_vertex_id ) ){
        return accessor->second;
    } else {
        ERROR("The given vertex does not exist: " << external_vertex_id);
    }
}

uint64_t LiveGraphDriver::int2ext(void* opaque_transaction, uint64_t internal_vertex_id) const {
    auto transaction = reinterpret_cast<lg::Transaction*>(opaque_transaction);
    string_view payload = transaction->get_vertex(internal_vertex_id);
    if(payload.empty()){ // the vertex does not exist
        return numeric_limits<uint64_t>::max();
    } else {
        return *(reinterpret_cast<const uint64_t*>(payload.data()));
    }
}

/*****************************************************************************
 *                                                                           *
 *  Updates                                                                  *
 *                                                                           *
 *****************************************************************************/
bool LiveGraphDriver::add_vertex(uint64_t external_id){
    //COUT_DEBUG("vertex_id: " << external_id);
    vertex_dictionary_t::accessor accessor; // xlock
    bool inserted = VertexDictionary->insert(accessor, external_id);
    if ( inserted ){
        lg::vertex_t internal_id = 0;
        bool done = false;
        do {
            try {
                auto tx = LiveGraph->begin_transaction();
                internal_id = tx.new_vertex();
                string_view data { (char*) &external_id, sizeof(external_id) };
                tx.put_vertex(internal_id, data);
                tx.commit();
                done = true;
            } catch(lg::Transaction::RollbackExcept& e){
                COUT_DEBUG("Rollback, vertex id: " << external_id);
                // retry ...
            }
        } while(!done);

        accessor->second = internal_id;
        m_num_vertices++;
    }

    return inserted;
}

bool LiveGraphDriver::remove_vertex(uint64_t external_id){
    COUT_DEBUG("vertex_id: " << external_id);
    vertex_dictionary_t::accessor accessor; // xlock
    bool found = VertexDictionary->find(accessor, external_id);
    if(found){
        lg::vertex_t internal_id = accessor->second;
        bool done = false;
        do {
            try {
                auto tx = LiveGraph->begin_transaction();
                tx.del_vertex(internal_id);
                tx.commit();
                done = true;
            } catch(lg::Transaction::RollbackExcept& e){
                // retry ...
            }
        } while(!done);

        VertexDictionary->erase(accessor);
    }
    m_num_vertices--;
    return found;
}

bool LiveGraphDriver::has_vertex(uint64_t vertex_id) const {
    vertex_dictionary_t::const_accessor accessor;
    return VertexDictionary->find(accessor, vertex_id);
}

bool LiveGraphDriver::add_edge(gfe::graph::WeightedEdge e){
    //COUT_DEBUG("Edge: " << e);

    vertex_dictionary_t::const_accessor accessor1, accessor2;  // shared lock on the dictionary
    if(!VertexDictionary->find(accessor1, e.source())){ return false; }
    if(!VertexDictionary->find(accessor2, e.destination())) { return false; }
    lg::vertex_t internal_source_id = accessor1->second;
    lg::vertex_t internal_destination_id = accessor2->second;

    bool done = false;
    do {
        try {
            auto tx = LiveGraph->begin_transaction();

            // insert the new edge only if it doesn't already exist
            auto lg_weight = tx.get_edge(internal_source_id, 0, internal_destination_id);
            if(lg_weight.size() > 0){ // the edge already exists
                tx.abort();
                return false;
            }

            string_view weight { (char*) &e.m_weight, sizeof(e.weight()) };
            tx.put_edge(internal_source_id, /* label */ 0, internal_destination_id, weight);
            if(!m_is_directed){ // undirected graph
                // We follow the same convention given by G. Feng, author of the LiveGraph paper,
                // for his experiments in the LDBC SNB Person knows Person: undirected edges
                // are added twice as a -> b and b -> a
                tx.put_edge(internal_destination_id, /* label */ 0, internal_source_id, weight);
            }

            tx.commit();
            m_num_edges++;
            done = true;
        } catch (lg::Transaction::RollbackExcept& exc){
            COUT_DEBUG("Rollback, edge: " << e);
            // retry ...
        }
    } while(!done);

    return true;
}

bool LiveGraphDriver::add_edge_v2(gfe::graph::WeightedEdge edge){
    uint64_t internal_source_id = numeric_limits<uint64_t>::max();
    uint64_t internal_destination_id = 0;
    bool insert_source = false;
    bool insert_destination = false;
    vertex_dictionary_t::const_accessor slock1, slock2;
    vertex_dictionary_t::accessor xlock1, xlock2;

    if(VertexDictionary->find(slock1, edge.m_source)){ // insert the vertex e.m_source
        internal_source_id = slock1->second;
    } else {
        slock1.release();
        if ( VertexDictionary->insert(xlock1, edge.m_source) ) {
            insert_source = true;
        } else {
            internal_source_id = xlock1->second;
        }
    }

    if(VertexDictionary->find(slock2, edge.m_destination)){ // insert the vertex e.m_destination
        internal_destination_id = slock2->second;
    } else {
        slock2.release();
        if( VertexDictionary->insert(xlock2, edge.m_destination) ) {
            insert_destination = true;
        } else {
            internal_destination_id = xlock2->second;
        }
    }

    bool done = false;
    do {
        try {
            auto tx = LiveGraph->begin_transaction();

            // create the vertices in LiveGraph
            if(insert_source){
                internal_source_id = tx.new_vertex();
                string_view data { (char*) &edge.m_source, sizeof(edge.m_source) };
                tx.put_vertex(internal_source_id, data);
            }
            if(insert_destination){
                internal_destination_id = tx.new_vertex();
                string_view data { (char*) &edge.m_destination, sizeof(edge.m_destination) };
                tx.put_vertex(internal_destination_id, data);
            }

            // insert the edge
            string_view weight { (char*) &edge.m_weight, sizeof(edge.m_weight) };
            tx.put_edge(internal_source_id, /* label */ 0, internal_destination_id, weight);

            // a) In directed graphs, we register the incoming edges with label 1
            // b) In undirected graphs, we follow the same convention given by G. Feng, author
            // of the LiveGraph paper, for his experiments in the LDBC SNB Person knows Person:
            // undirected edges are added twice as a -> b and b -> a
            lg::label_t label = m_is_directed ? 1 : 0;
            tx.put_edge(internal_destination_id, /* label */ label, internal_source_id, weight);

            tx.commit();
            m_num_edges++;

            done = true;
        } catch (lg::Transaction::RollbackExcept& e){
            // retry ...
        }
    } while(!done);

    if(insert_source){
        assert(internal_source_id != numeric_limits<uint64_t>::max());
        xlock1->second = internal_source_id;
        m_num_vertices++;
    }
    if(insert_destination){
        assert(internal_destination_id != numeric_limits<uint64_t>::max());
        xlock2->second = internal_destination_id;
       m_num_vertices++;
    }

    return true;
}

bool LiveGraphDriver::remove_edge(gfe::graph::Edge e){
    vertex_dictionary_t::const_accessor slock1, slock2;
    if(!VertexDictionary->find(slock1, e.source())){ return false; }
    if(!VertexDictionary->find(slock2, e.destination())){ return false; }
    lg::vertex_t internal_source_id = slock1->second;
    lg::vertex_t internal_destination_id = slock2->second;

    while(true){
        try {
            auto tx = LiveGraph->begin_transaction();
            bool removed = tx.del_edge(internal_source_id, /* label */ 0, internal_destination_id);
            if(removed && !m_is_directed){ // undirected graph
                tx.del_edge(internal_destination_id, /* label */ 0, internal_source_id);
            }
            tx.commit();
            if(removed){ m_num_edges--; }
            return removed;
        } catch(lg::Transaction::RollbackExcept& e){
            // retry ...
        }
    }
}

double LiveGraphDriver::get_weight(uint64_t source, uint64_t destination) const {
    // check whether the referred vertices exist
    vertex_dictionary_t::const_accessor slock1, slock2;
    if(!VertexDictionary->find(slock1, source)){ return numeric_limits<double>::signaling_NaN(); }
    if(!VertexDictionary->find(slock2, destination)){ return numeric_limits<double>::signaling_NaN(); }
    lg::vertex_t internal_source_id = slock1->second;
    lg::vertex_t internal_destination_id = slock2->second;

    auto tx = LiveGraph->begin_read_only_transaction();
    string_view lg_weight = tx.get_edge(internal_source_id, /* label */ 0, internal_destination_id);
    double weight = numeric_limits<double>::signaling_NaN();
    if(lg_weight.size() > 0){ // the edge exists
        weight = *(reinterpret_cast<const double*>(lg_weight.data()));
    }
    tx.abort(); // commit() fires the exception `The transaction is read-only without cache.'

    return weight;
}


/*****************************************************************************
 *                                                                           *
 *  Dump                                                                     *
 *                                                                           *
 *****************************************************************************/
void LiveGraphDriver::dump_ostream(std::ostream& out) const {
    out << "[LiveGraph] num vertices: " << m_num_vertices << ", num edges: " << m_num_edges << ", "
            "directed graph: " << boolalpha << is_directed() << ", read only txn for graphalytics: " << m_read_only << endl;
    auto tx = LiveGraph->begin_read_only_transaction();
    const uint64_t max_vertex_id = LiveGraph->get_max_vertex_id();
    for(uint64_t internal_source_id = 0; internal_source_id < max_vertex_id; internal_source_id++){
        auto lg_external_id = tx.get_vertex(internal_source_id);
        if(lg_external_id.size() == 0) continue; // the vertex has been deleted
        uint64_t external_id = *((uint64_t*) lg_external_id.data() );
        out << "[" << internal_source_id << ", external_id: " << external_id << "]";
        { // outgoing edges
            out << " outgoing edges: ";
            auto it = tx.get_edges(internal_source_id, 0);
            bool first = true;
            while(it.valid()){
                if(first){ first = false; } else { out << ", "; }
                uint64_t internal_destination_id = it.dst_id();
                auto lg_weight = it.edge_data();
                double weight = * ((double*) lg_weight.data());
                out << "<" << internal_destination_id << " [external: " << int2ext(&tx, internal_destination_id) << "], " << weight << ">";
                it.next();
            }
        }
        out << endl;
    }
    tx.abort(); // commit() fires the exception `The transaction is read-only without cache.'
}

/*****************************************************************************
 *                                                                           *
 *  Graphalytics Helpers                                                     *
 *                                                                           *
 *****************************************************************************/

template <typename T>
vector<pair<uint64_t, T>> LiveGraphDriver::translate(void* /* transaction object */ transaction, const T* __restrict data, uint64_t data_sz) {
    assert(transaction != nullptr && "Transaction object not specified");
    vector<pair<uint64_t, T>> output(data_sz);

    for(uint64_t logical_id = 0; logical_id < data_sz; logical_id++){
        uint64_t external_id = int2ext(transaction, logical_id);
        if(external_id == numeric_limits<uint64_t>::max()) { // the vertex does not exist
            output[logical_id] = make_pair(numeric_limits<uint64_t>::max(), numeric_limits<T>::max()); // special marker
        } else {
            output[logical_id] = make_pair(external_id, data[logical_id]);
        }
    }

    return output;
}

template <typename T, bool negative_scores>
void LiveGraphDriver::save_results(const vector<pair<uint64_t, T>>& result, const char* dump2file) {
    assert(dump2file != nullptr);
    COUT_DEBUG("save the results to: " << dump2file);

    fstream handle(dump2file, ios_base::out);
    if (!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

    for (const auto &p : result) {
        if(p.first == numeric_limits<uint64_t>::max()) continue; // invalid node

        handle << p.first << " ";

        if(!negative_scores && p.second < 0){
            handle << numeric_limits<T>::max();
        } else {
            handle << p.second;
        }

        handle << "\n";
    }

    handle.close();
}


/*****************************************************************************
 *                                                                           *
 *  BFS                                                                      *
 *                                                                           *
 *****************************************************************************/
// Implementation based on the reference BFS for the GAP Benchmark Suite
// https://github.com/sbeamer/gapbs
// The reference implementation has been written by Scott Beamer
//
// Copyright (c) 2015, The Regents of the University of California (Regents)
// All Rights Reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
// 3. Neither the name of the Regents nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL REGENTS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

/*

Will return parent array for a BFS traversal from a source vertex
This BFS implementation makes use of the Direction-Optimizing approach [1].
It uses the alpha and beta parameters to determine whether to switch search
directions. For representing the frontier, it uses a SlidingQueue for the
top-down approach and a Bitmap for the bottom-up approach. To reduce
false-sharing for the top-down approach, thread-local QueueBuffer's are used.
To save time computing the number of edges exiting the frontier, this
implementation precomputes the degrees in bulk at the beginning by storing
them in parent array as negative numbers. Thus the encoding of parent is:
  parent[x] < 0 implies x is unvisited and parent[x] = -out_degree(x)
  parent[x] >= 0 implies x been visited
[1] Scott Beamer, Krste Asanović, and David Patterson. "Direction-Optimizing
    Breadth-First Search." International Conference on High Performance
    Computing, Networking, Storage and Analysis (SC), Salt Lake City, Utah,
    November 2012.

*/
//#define DEBUG_BFS
#if defined(DEBUG_BFS)
#define COUT_DEBUG_BFS(msg) COUT_DEBUG(msg)
#else
#define COUT_DEBUG_BFS(msg)
#endif


static
int64_t do_bfs_BUStep(lg::Transaction& transaction, uint64_t max_vertex_id, int64_t* distances, int64_t distance, gapbs::Bitmap &front, gapbs::Bitmap &next) {
    int64_t awake_count = 0;
    next.reset();

    #pragma omp parallel for schedule(dynamic, 1024) reduction(+ : awake_count)
    for (uint64_t u = 0; u < max_vertex_id; u++) {
        if(distances[u] == numeric_limits<int64_t>::max()) continue; // the vertex does not exist
        COUT_DEBUG_BFS("explore: " << u << ", distance: " << distances[u]);

        if (distances[u] < 0){ // the node has not been visited yet
            auto iterator = transaction.get_edges(u, /* label */ 0); // fixme: incoming edges for directed graphs

            while(iterator.valid()){
                uint64_t dst = iterator.dst_id();
                COUT_DEBUG_BFS("\tincoming edge: " << dst);

                if(front.get_bit(dst)) {
                    COUT_DEBUG_BFS("\t-> distance updated to " << distance << " via vertex #" << dst);
                    distances[u] = distance; // on each BUStep, all nodes will have the same distance
                    awake_count++;
                    next.set_bit(u);
                    break;
                }

                iterator.next();
            }
        }
    }

    return awake_count;
}

static
int64_t do_bfs_TDStep(lg::Transaction& transaction, uint64_t max_vertex_id, int64_t* distances, int64_t distance, gapbs::SlidingQueue<int64_t>& queue) {
    int64_t scout_count = 0;

    #pragma omp parallel reduction(+ : scout_count)
    {
        gapbs::QueueBuffer<int64_t> lqueue(queue);

        #pragma omp for schedule(dynamic, 64)
        for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
            int64_t u = *q_iter;
            COUT_DEBUG_BFS("explore: " << u);
            auto iterator = transaction.get_edges(u, /* label */ 0);
            while(iterator.valid()){
                uint64_t dst = iterator.dst_id();
                COUT_DEBUG_BFS("\toutgoing edge: " << dst);

                int64_t curr_val = distances[dst];
                if (curr_val < 0 && gapbs::compare_and_swap(distances[dst], curr_val, distance)) {
                    COUT_DEBUG_BFS("\t-> distance updated to " << distance << " via vertex #" << dst);
                    lqueue.push_back(dst);
                    scout_count += -curr_val;
                }

                iterator.next();
            }
        }

        lqueue.flush();
    }

    return scout_count;
}

static
void do_bfs_QueueToBitmap(lg::Transaction& transaction, uint64_t max_vertex_id, const gapbs::SlidingQueue<int64_t> &queue, gapbs::Bitmap &bm) {
    #pragma omp parallel for
    for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
        int64_t u = *q_iter;
        bm.set_bit_atomic(u);
    }
}

static
void do_bfs_BitmapToQueue(lg::Transaction& transaction, uint64_t max_vertex_id, const gapbs::Bitmap &bm, gapbs::SlidingQueue<int64_t> &queue) {
    #pragma omp parallel
    {
        gapbs::QueueBuffer<int64_t> lqueue(queue);
        #pragma omp for
        for (uint64_t n=0; n < max_vertex_id; n++)
            if (bm.get_bit(n))
                lqueue.push_back(n);
        lqueue.flush();
    }
    queue.slide_window();
}

static
unique_ptr<int64_t[]> do_bfs_init_distances(lg::Transaction& transaction, uint64_t max_vertex_id) {
    unique_ptr<int64_t[]> distances{ new int64_t[max_vertex_id] };
    #pragma omp parallel for
    for (uint64_t n = 0; n < max_vertex_id; n++){
        if(transaction.get_vertex(n).empty()){ // the vertex does not exist
            distances[n] = numeric_limits<int64_t>::max();
        } else { // the vertex exists
            // Retrieve the out degree for the vertex n
            uint64_t out_degree = 0;
            auto iterator = transaction.get_edges(n, /* label */ 0);
            while(iterator.valid()){
                out_degree++;
                iterator.next();
            }

            distances[n] = out_degree != 0 ? - out_degree : -1;
        }
    }

    return distances;
}

static
unique_ptr<int64_t[]> do_bfs(lg::Transaction& transaction, uint64_t num_vertices, uint64_t num_edges, uint64_t max_vertex_id, uint64_t root, utility::TimeoutService& timer, int alpha = 15, int beta = 18) {
    // The implementation from GAP BS reports the parent (which indeed it should make more sense), while the one required by
    // Graphalytics only returns the distance
    unique_ptr<int64_t[]> ptr_distances = do_bfs_init_distances(transaction, max_vertex_id);
    int64_t* __restrict distances = ptr_distances.get();
    distances[root] = 0;

    gapbs::SlidingQueue<int64_t> queue(max_vertex_id);
    queue.push_back(root);
    queue.slide_window();
    gapbs::Bitmap curr(max_vertex_id);
    curr.reset();
    gapbs::Bitmap front(max_vertex_id);
    front.reset();
    int64_t edges_to_check = num_edges; //g.num_edges_directed();

    int64_t scout_count = 0;
    { // retrieve the out degree of the root
        auto iterator = transaction.get_edges(root, 0);
        while(iterator.valid()){ scout_count++; iterator.next(); }
    }
    int64_t distance = 1; // current distance

    while (!timer.is_timeout() && !queue.empty()) {

        if (scout_count > edges_to_check / alpha) {
            int64_t awake_count, old_awake_count;
            do_bfs_QueueToBitmap(transaction, max_vertex_id, queue, front);
            awake_count = queue.size();
            queue.slide_window();
            do {
                old_awake_count = awake_count;
                awake_count = do_bfs_BUStep(transaction, max_vertex_id, distances, distance, front, curr);
                front.swap(curr);
                distance++;
            } while ((awake_count >= old_awake_count) || (awake_count > (int64_t) num_vertices / beta));
            do_bfs_BitmapToQueue(transaction, max_vertex_id, front, queue);
            scout_count = 1;
        } else {
            edges_to_check -= scout_count;
            scout_count = do_bfs_TDStep(transaction, max_vertex_id, distances, distance, queue);
            queue.slide_window();
            distance++;
        }
    }

    return ptr_distances;
}

void LiveGraphDriver::bfs(uint64_t external_source_id, const char* dump2file) {
    if(m_is_directed) { ERROR("This implementation of the BFS does not support directed graphs"); }

    // Init
    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();
    lg::Transaction transaction = m_read_only ? LiveGraph->begin_read_only_transaction() : LiveGraph->begin_transaction();
    uint64_t max_vertex_id = LiveGraph->get_max_vertex_id();
    uint64_t num_vertices = m_num_vertices;
    uint64_t num_edges = m_num_edges;
    uint64_t root = ext2int(external_source_id);
    COUT_DEBUG_BFS("root: " << root << " [external vertex: " << external_source_id << "]");

    // Run the BFS algorithm
    unique_ptr<int64_t[]> ptr_result = do_bfs(transaction, num_vertices, num_edges, max_vertex_id, root, timeout);
    if(timeout.is_timeout()){
        transaction.abort(); // not sure if strictly necessary
        RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);
    }

    // translate the logical vertex IDs into the external vertex IDs
    auto external_ids = translate(&transaction, ptr_result.get(), max_vertex_id);
    transaction.abort(); // not sure if strictly necessary
    if(timeout.is_timeout()){
        RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);
    }

    if(dump2file != nullptr) // store the results in the given file
        save_results<int64_t, false>(external_ids, dump2file);
}

/*****************************************************************************
 *                                                                           *
 *  PageRank                                                                 *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG_PAGERANK
#if defined(DEBUG_PAGERANK)
#define COUT_DEBUG_PAGERANK(msg) COUT_DEBUG(msg)
#else
#define COUT_DEBUG_PAGERANK(msg)
#endif

// Implementation based on the reference PageRank for the GAP Benchmark Suite
// https://github.com/sbeamer/gapbs
// The reference implementation has been written by Scott Beamer
//
// Copyright (c) 2015, The Regents of the University of California (Regents)
// All Rights Reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
// 3. Neither the name of the Regents nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL REGENTS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

/*
GAP Benchmark Suite
Kernel: PageRank (PR)
Author: Scott Beamer

Will return pagerank scores for all vertices once total change < epsilon

This PR implementation uses the traditional iterative approach. This is done
to ease comparisons to other implementations (often use same algorithm), but
it is not necessarily the fastest way to implement it. It does perform the
updates in the pull direction to remove the need for atomics.
*/

static
unique_ptr<double[]> do_pagerank(lg::Transaction& transaction, uint64_t num_vertices, uint64_t max_vertex_id, uint64_t num_iterations, double damping_factor, utility::TimeoutService& timer) {
    const double init_score = 1.0 / num_vertices;
    const double base_score = (1.0 - damping_factor) / num_vertices;

    unique_ptr<double[]> ptr_scores{ new double[max_vertex_id]() }; // avoid memory leaks
    unique_ptr<uint64_t[]> ptr_degrees{ new uint64_t[max_vertex_id]() }; // avoid memory leaks
    double* scores = ptr_scores.get();
    uint64_t* __restrict degrees = ptr_degrees.get();


    #pragma omp parallel for
    for(uint64_t v = 0; v < max_vertex_id; v++){
        scores[v] = init_score;

        // compute the outdegree of the vertex
        if(!transaction.get_vertex(v).empty()){ // check the vertex exists
            uint64_t degree = 0;
            auto iterator = transaction.get_edges(v, /* label ? */ 0);
            while(iterator.valid()){ degree++; iterator.next(); }
            degrees[v] = degree;
        } else {
            degrees[v] = numeric_limits<uint64_t>::max();
        }
    }

    gapbs::pvector<double> outgoing_contrib(max_vertex_id, 0.0);

    // pagerank iterations
    for(uint64_t iteration = 0; iteration < num_iterations && !timer.is_timeout(); iteration++){
        double dangling_sum = 0.0;

        // for each node, precompute its contribution to all of its outgoing neighbours and, if it's a sink,
        // add its rank to the `dangling sum' (to be added to all nodes).
        #pragma omp parallel for reduction(+:dangling_sum)
        for(uint64_t v = 0; v < max_vertex_id; v++){
            uint64_t out_degree = degrees[v];
            if(out_degree == numeric_limits<uint64_t>::max()){
                continue; // the vertex does not exist
            } else if (out_degree == 0){ // this is a sink
                dangling_sum += scores[v];
            } else {
                outgoing_contrib[v] = scores[v] / out_degree;
            }
        }

        dangling_sum /= num_vertices;

        // compute the new score for each node in the graph
        #pragma omp parallel for schedule(dynamic, 64)
        for(uint64_t v = 0; v < max_vertex_id; v++){
            if(degrees[v] == numeric_limits<uint64_t>::max()){ continue; } // the vertex does not exist

            double incoming_total = 0;
            auto iterator = transaction.get_edges(v, /* label ? */ 0); // fixme: incoming edges for directed graphs
            while(iterator.valid()){
                uint64_t u = iterator.dst_id();
                incoming_total += outgoing_contrib[u];
                iterator.next();
            }

            // update the score
            scores[v] = base_score + damping_factor * (incoming_total + dangling_sum);
        }
    }

    return ptr_scores;
}

void LiveGraphDriver::pagerank(uint64_t num_iterations, double damping_factor, const char* dump2file) {
    if(m_is_directed) { ERROR("This implementation of PageRank does not support directed graphs"); }

    // Init
    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();
    lg::Transaction transaction = m_read_only? LiveGraph->begin_read_only_transaction() : LiveGraph->begin_transaction();
    uint64_t num_vertices = m_num_vertices;
    uint64_t max_vertex_id = LiveGraph->get_max_vertex_id();

    // Run the PageRank algorithm
    unique_ptr<double[]> ptr_result = do_pagerank(transaction, num_vertices, max_vertex_id, num_iterations, damping_factor, timeout);
    if(timeout.is_timeout()){ transaction.abort(); RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);  }

    // Retrieve the external node ids
    auto external_ids = translate(&transaction, ptr_result.get(), max_vertex_id);
    transaction.abort(); // read-only transaction, abort == commit
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // Store the results in the given file
    if(dump2file != nullptr)
        save_results(external_ids, dump2file);
}

/*****************************************************************************
 *                                                                           *
 *  WCC                                                                      *
 *                                                                           *
 *****************************************************************************/
// Implementation based on the reference WCC for the GAP Benchmark Suite
// https://github.com/sbeamer/gapbs
// The reference implementation has been written by Scott Beamer
//
// Copyright (c) 2015, The Regents of the University of California (Regents)
// All Rights Reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
// 3. Neither the name of the Regents nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL REGENTS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#define DEBUG_WCC
#if defined(DEBUG_WCC)
#define COUT_DEBUG_WCC(msg) COUT_DEBUG(msg)
#else
#define COUT_DEBUG_WCC(msg)
#endif

/*
GAP Benchmark Suite
Kernel: Connected Components (CC)
Author: Scott Beamer

Will return comp array labelling each vertex with a connected component ID

This CC implementation makes use of the Shiloach-Vishkin [2] algorithm with
implementation optimizations from Bader et al. [1]. Michael Sutton contributed
a fix for directed graphs using the min-max swap from [3], and it also produces
more consistent performance for undirected graphs.

[1] David A Bader, Guojing Cong, and John Feo. "On the architectural
    requirements for efficient execution of graph algorithms." International
    Conference on Parallel Processing, Jul 2005.

[2] Yossi Shiloach and Uzi Vishkin. "An o(logn) parallel connectivity algorithm"
    Journal of Algorithms, 3(1):57–67, 1982.

[3] Kishore Kothapalli, Jyothish Soman, and P. J. Narayanan. "Fast GPU
    algorithms for graph connectivity." Workshop on Large Scale Parallel
    Processing, 2010.
*/

// The hooking condition (comp_u < comp_v) may not coincide with the edge's
// direction, so we use a min-max swap such that lower component IDs propagate
// independent of the edge's direction.
static
unique_ptr<uint64_t[]> do_wcc(lg::Transaction& transaction, uint64_t max_vertex_id, utility::TimeoutService& timer) {
    // init
    COUT_DEBUG_WCC("max_vertex_id: " << max_vertex_id);
    unique_ptr<uint64_t[]> ptr_components { new uint64_t[max_vertex_id] };
    uint64_t* comp = ptr_components.get();

    #pragma omp parallel for
    for (uint64_t n = 0; n < max_vertex_id; n++){
        if(transaction.get_vertex(n).empty()){ // the vertex does not exist
            COUT_DEBUG_WCC("Vertex #" << n << " does not exist");
            comp[n] = numeric_limits<uint64_t>::max();
        } else {
            comp[n] = n;
        }
    }

    bool change = true;
    while (change && !timer.is_timeout()) {
        change = false;

        #pragma omp parallel for schedule(dynamic, 64)
        for (uint64_t u = 0; u < max_vertex_id; u++){
            if(comp[u] == numeric_limits<uint64_t>::max()) continue; // the vertex does not exist

            auto iterator = transaction.get_edges(u, 0);
            while(iterator.valid()){
                uint64_t v = iterator.dst_id();

                uint64_t comp_u = comp[u];
                uint64_t comp_v = comp[v];
                if (comp_u != comp_v) {
                    // Hooking condition so lower component ID wins independent of direction
                    uint64_t high_comp = std::max(comp_u, comp_v);
                    uint64_t low_comp = std::min(comp_u, comp_v);
                    if (high_comp == comp[high_comp]) {
                        change = true;
                        COUT_DEBUG_WCC("comp[" << high_comp << "] = " << low_comp);
                        comp[high_comp] = low_comp;
                    }
                }

                iterator.next();
            }
        }

        #pragma omp parallel for schedule(dynamic, 64)
        for (uint64_t n = 0; n < max_vertex_id; n++){
            if(comp[n] == numeric_limits<uint64_t>::max()) continue; // the vertex does not exist

            while (comp[n] != comp[comp[n]]) {
                comp[n] = comp[comp[n]];
            }
        }


        COUT_DEBUG_WCC("change: " << change);
    }

    return ptr_components;
}

void LiveGraphDriver::wcc(const char* dump2file) {
    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();
    auto transaction = m_read_only ? LiveGraph->begin_read_only_transaction() : LiveGraph->begin_transaction();
    uint64_t max_vertex_id = LiveGraph->get_max_vertex_id();

    // run wcc
    unique_ptr<uint64_t[]> ptr_components = do_wcc(transaction, max_vertex_id, timeout);
    if(timeout.is_timeout()){ transaction.abort(); RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // translate the vertex IDs
    auto external_ids = translate(&transaction, ptr_components.get(), max_vertex_id);
    transaction.abort(); // read-only transaction, abort == commit
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // store the results in the given file
    if(dump2file != nullptr)
        save_results(external_ids, dump2file);
}

/*****************************************************************************
 *                                                                           *
 *  CDLP                                                                     *
 *                                                                           *
 *****************************************************************************/
// same impl~ as the one done for llama
static
unique_ptr<uint64_t[]> do_cdlp(lg::Transaction& transaction, uint64_t max_vertex_id, bool is_graph_directed, uint64_t max_iterations, utility::TimeoutService& timer) {
   unique_ptr<uint64_t[]> ptr_labels0 { new uint64_t[max_vertex_id] };
   unique_ptr<uint64_t[]> ptr_labels1 { new uint64_t[max_vertex_id] };
   uint64_t* labels0 = ptr_labels0.get(); // current labels
   uint64_t* labels1 = ptr_labels1.get(); // labels for the next iteration

   // initialisation
   #pragma omp parallel for
   for(uint64_t v = 0; v < max_vertex_id; v++){
       string_view payload = transaction.get_vertex(v);
       if(payload.empty()){ // the vertex does not exist
           labels0[v] = labels1[v] = numeric_limits<uint64_t>::max();
       } else {
           labels0[v] = *reinterpret_cast<const uint64_t*>(payload.data());
       }
   }

   // algorithm pass
   bool change = true;
   uint64_t current_iteration = 0;
   while(current_iteration < max_iterations && change && !timer.is_timeout()){
       change = false; // reset the flag

       #pragma omp parallel for schedule(dynamic, 64) shared(change)
       for(uint64_t v = 0; v < max_vertex_id; v++){
           if(labels0[v] == numeric_limits<uint64_t>::max()) continue; // the vertex does not exist

           unordered_map<uint64_t, uint64_t> histogram;

           // compute the histogram from both the outgoing & incoming edges. The aim is to find the number of each label
           // is shared among the neighbours of node_id
           auto iterator = transaction.get_edges(v, 0); // out edges
           while(iterator.valid()){
               uint64_t u = iterator.dst_id();
               histogram[labels0[u]]++;
               iterator.next();
           }

           // get the max label
           uint64_t label_max = numeric_limits<int64_t>::max();
           uint64_t count_max = 0;
           for(const auto pair : histogram){
               if(pair.second > count_max || (pair.second == count_max && pair.first < label_max)){
                   label_max = pair.first;
                   count_max = pair.second;
               }
           }

           labels1[v] = label_max;
           change |= (labels0[v] != labels1[v]);
       }

       std::swap(labels0, labels1); // next iteration
       current_iteration++;
   }

   if(labels0 == ptr_labels0.get()){
       return ptr_labels0;
   } else {
       return ptr_labels1;
   }
}

void LiveGraphDriver::cdlp(uint64_t max_iterations, const char* dump2file) {
    if(m_is_directed) { ERROR("This implementation of the CDLP does not support directed graphs"); }

    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();
    auto transaction = m_read_only ?  LiveGraph->begin_read_only_transaction() : LiveGraph->begin_transaction();
    uint64_t max_vertex_id = LiveGraph->get_max_vertex_id();

    // Run the CDLP algorithm
    unique_ptr<uint64_t[]> labels = do_cdlp(transaction, max_vertex_id, is_directed(), max_iterations, timeout);
    if(timeout.is_timeout()){ transaction.abort(); RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);  }

    // Translate the vertex IDs
    auto external_ids = translate(&transaction, labels.get(), max_vertex_id);
    transaction.abort(); // read-only transaction, abort == commit
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // Store the results in the given file
    if(dump2file != nullptr)
        save_results(external_ids, dump2file);
}

/*****************************************************************************
 *                                                                           *
 *  LCC                                                                      *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG_LCC
#if defined(DEBUG_LCC)
#define COUT_DEBUG_LCC(msg) COUT_DEBUG(msg)
#else
#define COUT_DEBUG_LCC(msg)
#endif

// loosely based on the impl~ made for GraphOne
static
unique_ptr<double[]> do_lcc_undirected(lg::Transaction& transaction, uint64_t max_vertex_id, utility::TimeoutService& timer) {
    unique_ptr<double[]> ptr_lcc { new double[max_vertex_id] };
    double* lcc = ptr_lcc.get();
    unique_ptr<uint32_t[]> ptr_degrees_out { new uint32_t[max_vertex_id] };
    uint32_t* __restrict degrees_out = ptr_degrees_out.get();

    // precompute the degrees of the vertices
    #pragma omp parallel for schedule(dynamic, 4096)
    for(uint64_t v = 0; v < max_vertex_id; v++){
        bool vertex_exists = !transaction.get_vertex(v).empty();
        if(!vertex_exists){
            lcc[v] = numeric_limits<double>::signaling_NaN();
        } else {
            { // out degree, restrict the scope
                uint32_t count = 0;
                auto iterator = transaction.get_edges(v, 0);
                while(iterator.valid()){ count ++; iterator.next(); }
                degrees_out[v] = count;
            }
        }
    }


    #pragma omp parallel for schedule(dynamic, 64)
    for(uint64_t v = 0; v < max_vertex_id; v++){
        if(degrees_out[v] == numeric_limits<uint32_t>::max()) continue; // the vertex does not exist

        COUT_DEBUG_LCC("> Node " << v);
        if(timer.is_timeout()) continue; // exhausted the budget of available time
        lcc[v] = 0.0;
        uint64_t num_triangles = 0; // number of triangles found so far for the node v

        // Cfr. Spec v.0.9.0 pp. 15: "If the number of neighbors of a vertex is less than two, its coefficient is defined as zero"
        uint64_t v_degree_out = degrees_out[v];
        if(v_degree_out < 2) continue;

        // Build the list of neighbours of v
        unordered_set<uint64_t> neighbours;

        { // Fetch the list of neighbours of v
            auto iterator1 = transaction.get_edges(v, 0);
            while(iterator1.valid()){
                uint64_t u = iterator1.dst_id();
                neighbours.insert(u);
                iterator1.next();
            }
        }

        // again, visit all neighbours of v
        // for directed graphs, edges1 contains the intersection of both the incoming and the outgoing edges
        auto iterator1 = transaction.get_edges(v, /* label */ 0);
        while(iterator1.valid()){
            uint64_t u = iterator1.dst_id();
            COUT_DEBUG_LCC("[" << i << "/" << edges.size() << "] neighbour: " << u);
            assert(neighbours.count(u) == 1 && "The set `neighbours' should contain all neighbours of v");

            // For the Graphalytics spec v 0.9.0, only consider the outgoing edges for the neighbours u
            auto iterator2 = transaction.get_edges(u, /* label */ 0);
            while(iterator2.valid()){
                uint64_t w = iterator2.dst_id();

                COUT_DEBUG_LCC("---> [" << j << "/" << /* degree */ (u_out_interval.second - u_out_interval.first) << "] neighbour: " << w);
                // check whether it's also a neighbour of v
                if(neighbours.count(w) == 1){
                    COUT_DEBUG_LCC("Triangle found " << v << " - " << u << " - " << w);
                    num_triangles++;
                }

                iterator2.next();
            }

            iterator1.next();
        }

        // register the final score
        uint64_t max_num_edges = v_degree_out * (v_degree_out -1);
        lcc[v] = static_cast<double>(num_triangles) / max_num_edges;
        COUT_DEBUG_LCC("Score computed: " << (num_triangles) << "/" << max_num_edges << " = " << lcc[v]);
    }

    return ptr_lcc;
}

void LiveGraphDriver::lcc(const char* dump2file) {
    if(m_is_directed) { ERROR("Implementation of LCC supports only undirected graphs"); }

    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();
    lg::Transaction transaction = m_read_only ? LiveGraph->begin_read_only_transaction() : LiveGraph->begin_transaction();
    uint64_t max_vertex_id = LiveGraph->get_max_vertex_id();

    // Run the LCC algorithm
    unique_ptr<double[]> scores = do_lcc_undirected(transaction, max_vertex_id, timeout);
    if(timeout.is_timeout()){ transaction.abort(); RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);  }

    // Translate the vertex IDs
    auto external_ids = translate(&transaction, scores.get(), max_vertex_id);
    transaction.abort(); // read-only transaction, abort == commit
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // Store the results in the given file
    if(dump2file != nullptr)
        save_results(external_ids, dump2file);
}

/*****************************************************************************
 *                                                                           *
 *  SSSP                                                                     *
 *                                                                           *
 *****************************************************************************/
// Implementation based on the reference SSSP for the GAP Benchmark Suite
// https://github.com/sbeamer/gapbs
// The reference implementation has been written by Scott Beamer
//
// Copyright (c) 2015, The Regents of the University of California (Regents)
// All Rights Reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
// 3. Neither the name of the Regents nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL REGENTS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

using NodeID = uint64_t;
using WeightT = double;
static const size_t kMaxBin = numeric_limits<size_t>::max()/2;

static
gapbs::pvector<WeightT> do_sssp(lg::Transaction& transaction, uint64_t num_edges, uint64_t max_vertex_id, uint64_t source, double delta, utility::TimeoutService& timer) {
    // Init
    gapbs::pvector<WeightT> dist(max_vertex_id, numeric_limits<WeightT>::infinity());
    dist[source] = 0;
    gapbs::pvector<NodeID> frontier(num_edges);
    // two element arrays for double buffering curr=iter&1, next=(iter+1)&1
    size_t shared_indexes[2] = {0, kMaxBin};
    size_t frontier_tails[2] = {1, 0};
    frontier[0] = source;

    #pragma omp parallel
    {
        vector<vector<NodeID> > local_bins(0);
        size_t iter = 0;

        while (shared_indexes[iter&1] != kMaxBin) {
            size_t &curr_bin_index = shared_indexes[iter&1];
            size_t &next_bin_index = shared_indexes[(iter+1)&1];
            size_t &curr_frontier_tail = frontier_tails[iter&1];
            size_t &next_frontier_tail = frontier_tails[(iter+1)&1];
            #pragma omp for nowait schedule(dynamic, 64)
            for (size_t i=0; i < curr_frontier_tail; i++) {
                NodeID u = frontier[i];
                if (dist[u] >= delta * static_cast<WeightT>(curr_bin_index)) {
                    auto iterator = transaction.get_edges(u, /* label */ 0);
                    while(iterator.valid()){
                        uint64_t v = iterator.dst_id();
                        string_view payload = iterator.edge_data();
                        double w = *reinterpret_cast<const double*>(payload.data());

                        WeightT old_dist = dist[v];
                        WeightT new_dist = dist[u] + w;
                        if (new_dist < old_dist) {
                            bool changed_dist = true;
                            while (!gapbs::compare_and_swap(dist[v], old_dist, new_dist)) {
                                old_dist = dist[v];
                                if (old_dist <= new_dist) {
                                    changed_dist = false;
                                    break;
                                }
                            }
                            if (changed_dist) {
                                size_t dest_bin = new_dist/delta;
                                if (dest_bin >= local_bins.size()) {
                                    local_bins.resize(dest_bin+1);
                                }
                                local_bins[dest_bin].push_back(v);
                            }
                        }

                        iterator.next();
                    }
                }
            }

            for (size_t i=curr_bin_index; i < local_bins.size(); i++) {
                if (!local_bins[i].empty()) {
                    #pragma omp critical
                    next_bin_index = min(next_bin_index, i);
                    break;
                }
            }

            #pragma omp barrier
            #pragma omp single nowait
            {
                curr_bin_index = kMaxBin;
                curr_frontier_tail = 0;
            }

            if (next_bin_index < local_bins.size()) {
                size_t copy_start = gapbs::fetch_and_add(next_frontier_tail, local_bins[next_bin_index].size());
                copy(local_bins[next_bin_index].begin(), local_bins[next_bin_index].end(), frontier.data() + copy_start);
                local_bins[next_bin_index].resize(0);
            }

            iter++;
            #pragma omp barrier
        }

#if defined(DEBUG)
        #pragma omp single
        COUT_DEBUG("took " << iter << " iterations");
#endif
    }

    return dist;
}

void LiveGraphDriver::sssp(uint64_t source_vertex_id, const char* dump2file) {
    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();
    lg::Transaction transaction = m_read_only ? LiveGraph->begin_read_only_transaction() : LiveGraph->begin_transaction();
    uint64_t num_edges = m_num_edges;
    uint64_t max_vertex_id = LiveGraph->get_max_vertex_id();
    uint64_t root = ext2int(source_vertex_id);

    // Run the SSSP algorithm
    double delta = 2.0; // same value used in the GAPBS, at least for most graphs
    auto distances = do_sssp(transaction, num_edges, max_vertex_id, root, delta, timeout);
    if(timeout.is_timeout()){ transaction.abort(); RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);  }

    // Translate the vertex IDs
    auto external_ids = translate(&transaction, distances.data(), max_vertex_id);
    transaction.abort(); // read-only transaction, abort == commit
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // Store the results in the given file
    if(dump2file != nullptr)
        save_results(external_ids, dump2file);
}

} // namespace
