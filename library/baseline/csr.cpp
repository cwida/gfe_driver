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

#include "csr.hpp"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <mutex>
#if defined(HAVE_LIBNUMA)
#include <numa.h>
#endif
#include <omp.h>
#include <random>
#include <sstream>
#include <string>
#include <unordered_set>

#include "common/error.hpp"
#include "common/system.hpp"
#include "common/timer.hpp"
#include "graph/edge_stream.hpp"
#include "graph/vertex_list.hpp"
#include "third-party/gapbs/gapbs.hpp"
#include "third-party/libcuckoo/cuckoohash_map.hh"
#include "utility/timeout_service.hpp"

using namespace common;
using namespace libcuckoo;
using namespace std;

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG
namespace gfe { extern mutex _log_mutex [[maybe_unused]]; }
#define COUT_DEBUG_FORCE(msg) { std::scoped_lock<std::mutex> lock{::gfe::_log_mutex}; std::cout << "[CSR::" << __FUNCTION__ << "] [Thread #" << common::concurrency::get_thread_id() << "] " << msg << std::endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif

/*****************************************************************************
 *                                                                           *
 *  Initialisation                                                           *
 *                                                                           *
 *****************************************************************************/

namespace gfe::library {

CSR::CSR(bool is_directed, bool numa_interleaved) : m_is_directed(is_directed), m_num_vertices (0), m_num_edges(0), m_numa_interleaved(numa_interleaved) {
#if !defined(HAVE_LIBNUMA)
    ERROR("[CSR] Cannot allocate the memory interleaved, dependency on libnuma missing");
#else
    if(numa_available() < 0){
        ERROR("[CSR] Cannot allocate the memory interleaved, a call to numa_available() returns a negative value (=> NUMA not available)");
    }
#endif
}

CSR::~CSR(){
    free_array(m_out_v); m_out_v = nullptr;
    free_array(m_out_e); m_out_e = nullptr;
    free_array(m_out_w); m_out_w = nullptr;

    if(m_is_directed){ // otherwise, they are simply aliases to m_out_x
        free_array(m_in_v);
        free_array(m_in_e);
        free_array(m_in_w);
    }

    m_in_v = nullptr;
    m_in_e = nullptr;
    m_in_w = nullptr;

    free_array(m_log2ext); m_log2ext = nullptr;
}

template<typename T>
T* CSR::alloca_array(uint64_t array_sz){
    if(m_numa_interleaved){
#if defined(HAVE_LIBNUMA)
        uint64_t required_bytes = /* header */ sizeof(uint64_t) + /* data */ sizeof(T) * array_sz;
        void* ptr = numa_alloc_interleaved(required_bytes);
        if(ptr == nullptr){ ERROR("[CSR] Cannot allocate an interleaved array of " << array_sz * sizeof(T) << " bytes"); }
        memset(ptr, '\0', required_bytes);
        uint64_t* header = reinterpret_cast<uint64_t*>(ptr);
        header[0] = required_bytes;
        return reinterpret_cast<T*>(header + 1);
#else
        return nullptr;
#endif
    } else {
        return new T[array_sz]();
    }
}

template<typename T>
void CSR::free_array(T* array){
    if(array == nullptr) return; // nop

    if(m_numa_interleaved){
#if defined(HAVE_LIBNUMA)
        uint64_t* start = reinterpret_cast<uint64_t*>(array) -1;
        uint64_t allocation_size = start[0];
        numa_free(start, allocation_size);
#else
        return; // nop
#endif
    } else {
        delete[] array;
    }
}

/*****************************************************************************
 *                                                                           *
 *  Properties                                                               *
 *                                                                           *
 *****************************************************************************/

uint64_t CSR::num_edges() const {
    return m_num_edges;
}

uint64_t CSR::num_vertices() const {
    return m_num_vertices;
}

bool CSR::is_directed() const {
    return m_is_directed;
}

bool CSR::has_vertex(uint64_t vertex_id) const {
    return m_ext2log.count(vertex_id);
}

double CSR::get_weight(uint64_t source, uint64_t destination) const {
    uint64_t logical_source_id = 0, logical_destination_id = 0;
    try {
        logical_source_id = m_ext2log.at(source);
        logical_destination_id = m_ext2log.at(destination);
    } catch(out_of_range& e){ // either source or destination do not exist
        return numeric_limits<double>::signaling_NaN();
    }

    auto offset = get_out_interval(logical_source_id);
    for(uint64_t i = offset.first, end = offset.second; i < end && m_out_e[i] <= logical_destination_id; i++){
        if(m_out_e[i] == logical_destination_id){
            return m_out_w[i];
        }
    }

    return numeric_limits<double>::signaling_NaN();
}

pair<uint64_t, uint64_t> CSR::get_out_interval(uint64_t logical_vertex_id) const {
    return get_interval_impl(m_out_v, logical_vertex_id);
}

pair<uint64_t, uint64_t> CSR::get_in_interval(uint64_t logical_vertex_id) const {
    return get_interval_impl(m_in_v, logical_vertex_id);
}

pair<uint64_t, uint64_t> CSR::get_interval_impl(const uint64_t* __restrict vertex_array, uint64_t logical_vertex_id) const {
    assert(logical_vertex_id < m_num_vertices && "Invalid vertex ID");
    if(logical_vertex_id == 0){
        return make_pair(0ull, vertex_array[0]);
    } else {
        return make_pair(vertex_array[logical_vertex_id -1], vertex_array[logical_vertex_id]);
    }
}

uint64_t CSR::get_out_degree(uint64_t logical_vertex_id) const {
    auto interval = get_out_interval(logical_vertex_id);
    return interval.second - interval.first;
}

uint64_t CSR::get_in_degree(uint64_t logical_vertex_id) const {
    auto interval = get_in_interval(logical_vertex_id);
    return interval.second - interval.first;
}

uint64_t CSR::get_random_vertex_id() const {
    std::mt19937_64 generator { /* seed */ std::random_device{}() };
    std::uniform_int_distribution<uint64_t> distribution{ 0, m_num_vertices -1 };
    uint64_t outcome = distribution(generator);
    return m_log2ext[outcome];
}

void CSR::set_timeout(uint64_t seconds) {
    m_timeout = seconds;
}

uint64_t* CSR::out_v() const { return m_out_v; }
uint64_t* CSR::out_e() const { return m_out_e; }
double* CSR::out_w() const { return m_out_w; }
uint64_t* CSR::in_v() const { return m_in_v; }
uint64_t* CSR::in_e() const { return m_in_e; }
double* CSR::in_w() const { return m_in_w; }

/*****************************************************************************
 *                                                                           *
 *  Load                                                                     *
 *                                                                           *
 *****************************************************************************/
void CSR::load(const std::string& path){
    if(m_out_v != nullptr) ERROR("Already initialised & loaded");

    ::gfe::graph::WeightedEdgeStream stream { path };
    load(stream);
}

void CSR::load(gfe::graph::WeightedEdgeStream& stream){
    if(m_out_v != nullptr) ERROR("Already initialised & loaded");

    if(m_is_directed){
        load_directed(stream);
    } else {
        load_undirected(stream);
    }
}

void CSR::load_directed(gfe::graph::WeightedEdgeStream& stream){
    m_num_edges = stream.num_edges();

    { // init the mapping of vertices
        auto vertices = stream.vertex_list();
        vertices->sort();
        m_num_vertices = vertices->num_vertices();

        m_ext2log.reserve(m_num_vertices);
        m_log2ext = alloca_array<uint64_t>(m_num_vertices);

        for(uint64_t i = 0; i < m_num_vertices; i++){
            uint64_t vertex_id = vertices->get(i);
            m_ext2log[vertex_id] = i;
            m_log2ext[i] = vertex_id;
        }
    }

    // init the outgoing edges
    stream.sort_by_src_dst();
    m_out_v = alloca_array<uint64_t>(m_num_vertices); // init to 0
    m_out_e = alloca_array<uint64_t>(m_num_edges);
    m_out_w = alloca_array<double>(m_num_edges);

    { // restrict the scope
        uint64_t external_source_id = numeric_limits<uint64_t>::max();
        uint64_t logical_source_id = 0;
        for(uint64_t i = 0; i < m_num_edges; i++){
            auto edge = stream.get(i);
            if(edge.source() != external_source_id){
                external_source_id = edge.source();
                assert(m_ext2log.count(external_source_id) == 1 && "The source vertex is not registered in the mapping");
                logical_source_id = m_ext2log[external_source_id];
            }
            m_out_v[logical_source_id]++;

            assert(m_ext2log.count(edge.destination()) == 1 && "The destination vertex is not registered in the mapping");
            uint64_t logical_destination_id = m_ext2log[edge.destination()];

            m_out_e[i] = logical_destination_id;
            m_out_w[i] = edge.weight();
        }
    }

    // prefix sum on the vertex array
    for(uint64_t i = 1; i < m_num_vertices; i++){
        m_out_v[i] += m_out_v[i -1];
    }

    // init the incoming edges
    // do it for both directed and undirected graphs. For undirected graphs it's an intermediate step
    stream.sort_by_dst_src();
    m_in_v = alloca_array<uint64_t>(m_num_vertices); // init to 0
    m_in_e = alloca_array<uint64_t>(m_num_edges);
    m_in_w = alloca_array<double>(m_num_edges);

    uint64_t external_destination_id = numeric_limits<uint64_t>::max();
    uint64_t logical_destination_id = 0;

    for(uint64_t i = 0; i < m_num_edges; i++){
        auto edge = stream.get(i);
        if(edge.destination() != external_destination_id){
            external_destination_id = edge.destination();
            assert(m_ext2log.count(external_destination_id) == 1 && "The destination vertex is not registered in the mapping");
            logical_destination_id = m_ext2log[external_destination_id];
        }
        m_in_v[logical_destination_id]++;

        assert(m_ext2log.count(edge.source()) == 1 && "The source vertex is not registered in the mapping");
        uint64_t logical_source_id = m_ext2log[edge.source()];

        m_in_e[i] = logical_source_id;
        m_in_w[i] = edge.weight();
    }

    // prefix sum on the vertex array
    for(uint64_t i = 1; i < m_num_vertices; i++){
        m_in_v[i] += m_in_v[i -1];
    }
}

void CSR::load_undirected(gfe::graph::WeightedEdgeStream& stream){
    // We rely to load_directed to build a directed graph first
    load_directed(stream);

    // init the final vectors
    auto fn_free_array = [this](uint64_t* ptr){ this->free_array(ptr); };
    unique_ptr<uint64_t, decltype(fn_free_array)> ptr_out_v_undirected { alloca_array<uint64_t>(m_num_vertices), fn_free_array };
    unique_ptr<uint64_t, decltype(fn_free_array)> ptr_out_e_undirected { alloca_array<uint64_t>(m_num_edges *2), fn_free_array };
    auto fn_free_array_dbl = [this](double* ptr){ this->free_array(ptr); };
    unique_ptr<double, decltype(fn_free_array_dbl)> ptr_out_w_undirected { alloca_array<double>(m_num_edges *2), fn_free_array_dbl };
    auto out_v_undirected = ptr_out_v_undirected.get();
    auto out_e_undirected = ptr_out_e_undirected.get();
    auto out_w_undirected = ptr_out_w_undirected.get();

    // process all vertices from the directed CSR
    uint64_t k = 0; // index on the edge vector
    for(uint64_t v = 0; v < m_num_vertices; v++){
        auto interval_out = get_out_interval(v);
        auto interval_in = get_in_interval(v);
        auto degree_out = interval_out.second - interval_out.first;
        auto degree_in = interval_in.second - interval_in.first;
        out_v_undirected[v] = degree_out + degree_in;

        // merge the incoming & outcoming edges for v
        uint64_t out = interval_out.first;
        uint64_t in = interval_in.first;
        while(out < interval_out.second && in < interval_in.second){
            assert(m_out_e[out] != m_in_e[in] && "We don't expect the same edge in both directions");

            if(m_out_e[out] < m_in_e[in]){
                out_e_undirected[k] = m_out_e[out];
                out_w_undirected[k] = m_out_w[out];

                out++;
                k++;
            } else {
                out_e_undirected[k] = m_in_e[in];
                out_w_undirected[k] = m_in_w[in];

                in++;
                k++;
            }
        }
        while(out < interval_out.second){
            out_e_undirected[k] = m_out_e[out];
            out_w_undirected[k] = m_out_w[out];

            out++;
            k++;
        }
        while(in < interval_in.second){
            out_e_undirected[k] = m_in_e[in];
            out_w_undirected[k] = m_in_w[in];

            in++;
            k++;
        }
    }

    // prefix sum
    for(uint64_t v = 1; v < m_num_vertices; v++){
        out_v_undirected[v] += out_v_undirected[v -1];
    }

    // replace the out & in vectors
    free_array(m_out_v);
    free_array(m_out_e);
    free_array(m_out_w);
    free_array(m_in_v);
    free_array(m_in_e);
    free_array(m_in_w);
    m_out_v = m_in_v = ptr_out_v_undirected.release();
    m_out_e = m_in_e = ptr_out_e_undirected.release();
    m_out_w = m_in_w = ptr_out_w_undirected.release();
}

/*****************************************************************************
 *                                                                           *
 *  Dump                                                                     *
 *                                                                           *
 *****************************************************************************/

void CSR::dump_ostream(std::ostream& out) const {
    out << "[CSR] directed graph: " << (m_is_directed ? "yes" : "no");
    out << ", num vertices: " << m_num_vertices << ", num edges: " << m_num_edges << "\n";
    for(uint64_t logical_vertex_id = 0; logical_vertex_id < m_num_vertices; logical_vertex_id ++ ){
        stringstream ss;
        ss << "[" << logical_vertex_id << "] vtx: " << m_log2ext[logical_vertex_id] << ", ";
        uint64_t ss_len = ss.tellp();
        out << ss.str();
        out << "out edges: ";
        auto offset = get_out_interval(logical_vertex_id);
        auto degree = offset.second - offset.first;
        for(uint64_t j = 0; j < degree; j++){
            if(j > 0) out << ", ";
            out << "<" << m_log2ext[ m_out_e[offset.first + j] ] << " (logical: " << m_out_e[offset.first + j] << "), " << m_out_w[offset.first + j] << ">";
        }
        out << "\n";
        if(m_is_directed){
            printf("%*c", (int) ss_len, ' ');
            out << "in edges: ";
            offset = get_in_interval(logical_vertex_id);
            degree = offset.second - offset.first;
            for(uint64_t j = 0; j < degree; j++){
                if(j > 0) out << ", ";
                out << "<" << m_log2ext[ m_in_e[offset.first + j] ] << " (logical: " << m_in_e[offset.first + j] << "), " << m_in_w[offset.first + j] << ">";
            }
            out << "\n";
        }
    }
}

/*****************************************************************************
 *                                                                           *
 *  Utility                                                                  *
 *                                                                           *
 *****************************************************************************/
template <typename T>
vector<pair<uint64_t, T>> CSR::translate(const T* __restrict values, uint64_t N) {
    vector<pair<uint64_t , T>> logical_result(N);

    #pragma omp parallel for
    for (uint64_t v = 0; v < N; v++) {
        logical_result[v] = make_pair(m_log2ext[v], values[v]);
    }
    return logical_result;
}

template <typename T, bool negative_scores>
void CSR::save_results(const vector<pair<uint64_t, T>>& result, const char* dump2file) {
    assert(dump2file != nullptr);
    COUT_DEBUG("save the results to: " << dump2file);

    fstream handle(dump2file, ios_base::out);
    if (!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

    for (const auto &p : result) {
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


int64_t CSR::do_bfs_BUStep(int64_t* distances, int64_t distance, gapbs::Bitmap &front, gapbs::Bitmap &next) const {
    int64_t awake_count = 0;
    next.reset();
    uint64_t* __restrict in_e = m_in_e;

    #pragma omp parallel for schedule(dynamic, 1024) reduction(+ : awake_count)
    for (uint64_t u = 0; u < m_num_vertices; u++) {
        COUT_DEBUG_BFS("explore " << u << " [external vertex = " << m_log2ext[u] << "], distance: " << distances[u]);

        if (distances[u] < 0){ // the node has not been visited yet
            auto in_interval = get_in_interval(u);
            uint64_t degree = in_interval.second - in_interval.first;
            if(degree == 0) continue;

            for(uint64_t i = in_interval.first; i < in_interval.second; i++){
                uint64_t dst = in_e[i];
                COUT_DEBUG_BFS("\tincoming edge: " << dst << " [external vertex = " << m_log2ext[dst] << "]");

                if(front.get_bit(dst)) {
                    COUT_DEBUG_BFS("\t-> distance updated to " << distance << " via vertex #" << dst << " [external vertex = " << m_log2ext[dst] << "]");
                    distances[u] = distance; // on each BUStep, all nodes will have the same distance
                    awake_count++;
                    next.set_bit(u);
                    break;
                }
            }
        }
    }

    return awake_count;
}

int64_t CSR::do_bfs_TDStep(int64_t* distances, int64_t distance, gapbs::SlidingQueue<int64_t>& queue) const {
    int64_t scout_count = 0;
    uint64_t* __restrict out_e = m_out_e;

    #pragma omp parallel reduction(+ : scout_count)
    {
        gapbs::QueueBuffer<int64_t> lqueue(queue);

        #pragma omp for schedule(dynamic, 64)
        for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
            int64_t u = *q_iter;
            COUT_DEBUG_BFS("explore: " << u << " [external vertex = "  << m_log2ext[u] << "]");
            auto out_interval = get_out_interval(u);
            uint64_t degree = out_interval.second - out_interval.first;
            if(degree == 0) continue;

            for(uint64_t i = out_interval.first; i < out_interval.second; i++){
                uint64_t dst = out_e[i];
                COUT_DEBUG_BFS("\toutgoing edge: " << dst << " [external vertex = " << m_log2ext[dst] << "]");

                int64_t curr_val = distances[dst];
                if (curr_val < 0 && gapbs::compare_and_swap(distances[dst], curr_val, distance)) {
                    COUT_DEBUG_BFS("\t-> distance updated to " << distance << " via vertex #" << dst << " [external vertex = " << m_log2ext[dst] << "]");
                    lqueue.push_back(dst);
                    scout_count += -curr_val;
                }
            }
        }

        lqueue.flush();
    }

    return scout_count;
}

void CSR::do_bfs_QueueToBitmap(const gapbs::SlidingQueue<int64_t> &queue, gapbs::Bitmap &bm) const {
    #pragma omp parallel for
    for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
        int64_t u = *q_iter;
        bm.set_bit_atomic(u);
    }
}

void CSR::do_bfs_BitmapToQueue(const gapbs::Bitmap &bm, gapbs::SlidingQueue<int64_t> &queue) const {
    #pragma omp parallel
    {
        gapbs::QueueBuffer<int64_t> lqueue(queue);
        #pragma omp for
        for (uint64_t n=0; n < m_num_vertices; n++)
            if (bm.get_bit(n))
                lqueue.push_back(n);
        lqueue.flush();
    }
    queue.slide_window();
}

unique_ptr<int64_t[]> CSR::do_bfs_init_distances() const {
    unique_ptr<int64_t[]> distances{ new int64_t[m_num_vertices] };
    #pragma omp parallel for
    for (uint64_t n = 0; n < m_num_vertices; n++){
        int64_t out_degree = get_out_degree(n);
        distances[n] = out_degree != 0 ? - out_degree : -1;
    }
    return distances;
}

unique_ptr<int64_t[]> CSR::do_bfs(uint64_t root, utility::TimeoutService& timer, int alpha, int beta) const {
    // The implementation from GAP BS reports the parent (which indeed it should make more sense), while the one required by
    // Graphalytics only returns the distance
    unique_ptr<int64_t[]> ptr_distances = do_bfs_init_distances();
    int64_t* __restrict distances = ptr_distances.get();
    distances[root] = 0;

    gapbs::SlidingQueue<int64_t> queue(m_num_vertices);
    queue.push_back(root);
    queue.slide_window();
    gapbs::Bitmap curr(m_num_vertices);
    curr.reset();
    gapbs::Bitmap front(m_num_vertices);
    front.reset();
    int64_t edges_to_check = m_num_edges; //g.num_edges_directed();
    int64_t scout_count = get_out_degree(root);
    int64_t distance = 1; // current distance
    while (!timer.is_timeout() && !queue.empty()) {

        if (scout_count > edges_to_check / alpha) {
            int64_t awake_count, old_awake_count;
            do_bfs_QueueToBitmap(queue, front);
            awake_count = queue.size();
            queue.slide_window();
            do {
                old_awake_count = awake_count;
                awake_count = do_bfs_BUStep(distances, distance, front, curr);
                front.swap(curr);
                distance++;
            } while ((awake_count >= old_awake_count) || (awake_count > (int64_t) m_num_vertices / beta));
            do_bfs_BitmapToQueue(front, queue);
            scout_count = 1;
        } else {
            edges_to_check -= scout_count;
            scout_count = do_bfs_TDStep(distances, distance, queue);
            queue.slide_window();
            distance++;
        }
    }

    return ptr_distances;
}

void CSR::bfs(uint64_t external_source_id, const char* dump2file) {
    // Init
    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();
    uint64_t root = m_ext2log.at(external_source_id);
    COUT_DEBUG_BFS("root: " << root << " [external vertex: " << external_source_id << "]");

    // Run the BFS algorithm
    unique_ptr<int64_t[]> ptr_result = do_bfs(root, timeout);
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);  }

    // Translate the logical IDs into the external IDs
    auto translation = translate(ptr_result.get(), m_num_vertices);
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);  }

    // Store the results in the given file
    if(dump2file != nullptr)
        save_results<int64_t, false>(translation, dump2file);
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
it is not necesarily the fastest way to implement it. It does perform the
updates in the pull direction to remove the need for atomics.
*/

unique_ptr<double[]> CSR::do_pagerank(uint64_t num_iterations, double damping_factor, utility::TimeoutService& timer) const {
    const double init_score = 1.0 / m_num_vertices;
    const double base_score = (1.0 - damping_factor) / m_num_vertices;

    unique_ptr<double[]> ptr_scores{ new double[m_num_vertices]() }; // avoid memory leaks
    double* scores = ptr_scores.get();
    #pragma omp parallel for
    for(uint64_t v = 0; v < m_num_vertices; v++){
        scores[v] = init_score;
    }
    gapbs::pvector<double> outgoing_contrib(m_num_vertices, 0.0);
    uint64_t* __restrict in_e = m_in_e;

    // pagerank iterations
    for(uint64_t iteration = 0; iteration < num_iterations && !timer.is_timeout(); iteration++){
        double dangling_sum = 0.0;

        // for each node, precompute its contribution to all of its outgoing neighbours and, if it's a sink,
        // add its rank to the `dangling sum' (to be added to all nodes).
        #pragma omp parallel for reduction(+:dangling_sum)
        for(uint64_t v = 0; v < m_num_vertices; v++){
            uint64_t out_degree = get_out_degree(v);
            if(out_degree == 0){ // this is a sink
                dangling_sum += scores[v];
            } else {
                outgoing_contrib[v] = scores[v] / out_degree;
            }
        }

        dangling_sum /= m_num_vertices;

        // compute the new score for each node in the graph
        #pragma omp parallel for schedule(dynamic, 64)
        for(uint64_t v = 0; v < m_num_vertices; v++){
            auto in_interval = get_in_interval(v);
            double incoming_total = 0;
            for(uint64_t i = in_interval.first; i < in_interval.second; i++){
                uint64_t u = in_e[i];
                incoming_total += outgoing_contrib[u];
            }

            // update the score
            scores[v] = base_score + damping_factor * (incoming_total + dangling_sum);
        }
    }

    return ptr_scores;
}

void CSR::pagerank(uint64_t num_iterations, double damping_factor, const char* dump2file) {
    // Init
    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();

    // Run the PageRank algorithm
    unique_ptr<double[]> ptr_result = do_pagerank(num_iterations, damping_factor, timeout);
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);  }

    // retrieve the external node ids
    auto translation = translate(ptr_result.get(), m_num_vertices);
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // store the results in the given file
    if(dump2file != nullptr){
        save_results(translation, dump2file);
    }
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


//#define DEBUG_WCC
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
unique_ptr<uint64_t[]> CSR::do_wcc(utility::TimeoutService& timer) const {
    // init
    unique_ptr<uint64_t[]> ptr_components { new uint64_t[m_num_vertices] };
    uint64_t* comp = ptr_components.get();
    uint64_t* __restrict out_e = m_out_e;

    #pragma omp parallel for
    for (uint64_t n = 0; n < m_num_vertices; n++){
        comp[n] = n;
    }

    bool change = true;
    while (change && !timer.is_timeout()) {
        change = false;

        #pragma omp parallel for schedule(dynamic, 64)
        for (uint64_t u = 0; u < m_num_vertices; u++){
            auto out_interval = get_out_interval(u);
            for(uint64_t i = out_interval.first; i < out_interval.second; i++){
                uint64_t v = out_e[i];

                uint64_t comp_u = comp[u];
                uint64_t comp_v = comp[v];
                if (comp_u == comp_v) continue;
                // Hooking condition so lower component ID wins independent of direction
                uint64_t high_comp = std::max(comp_u, comp_v);
                uint64_t low_comp = std::min(comp_u, comp_v);
                if (high_comp == comp[high_comp]) {
                    change = true;
                    comp[high_comp] = low_comp;
                }
            }
        }

        #pragma omp parallel for schedule(dynamic, 64)
        for (uint64_t n = 0; n < m_num_vertices; n++){
            while (comp[n] != comp[comp[n]]) {
                comp[n] = comp[comp[n]];
            }
        }
    }

    return ptr_components;
}

void CSR::wcc(const char* dump2file) {
    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();

    // run wcc
    unique_ptr<uint64_t[]> ptr_components = do_wcc(timeout);

    // retrieve the external node ids
    auto translation = translate(ptr_components.get(), m_num_vertices);
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // store the results in the given file
    if(dump2file != nullptr){
        save_results(translation, dump2file);
    }
}

/*****************************************************************************
 *                                                                           *
 *  CDLP                                                                     *
 *                                                                           *
 *****************************************************************************/
// same impl~ as the one done for llama
unique_ptr<uint64_t[]> CSR::do_cdlp(uint64_t max_iterations, utility::TimeoutService& timer) const {
   unique_ptr<uint64_t[]> ptr_labels0 { new uint64_t[m_num_vertices] };
   unique_ptr<uint64_t[]> ptr_labels1 { new uint64_t[m_num_vertices] };
   uint64_t* labels0 = ptr_labels0.get(); // current labels
   uint64_t* labels1 = ptr_labels1.get(); // labels for the next iteration

   // initialisation
   #pragma omp parallel for
   for(uint64_t v = 0; v < m_num_vertices; v++){
       labels0[v] = m_log2ext[v];
   }

   // algorithm pass
   bool change = true;
   uint64_t current_iteration = 0;
   uint64_t* __restrict out_e = m_out_e;
   uint64_t* __restrict in_e = m_in_e;
   while(current_iteration < max_iterations && change && !timer.is_timeout()){
       change = false; // reset the flag

       #pragma omp parallel for schedule(dynamic, 64) shared(change)
       for(uint64_t v = 0; v < m_num_vertices; v++){
           unordered_map<uint64_t, uint64_t> histogram;

           // compute the histogram from both the outgoing & incoming edges. The aim is to find the number of each label
           // is shared among the neighbours of node_id
           auto out_interval = get_out_interval(v);
           for(uint64_t i = out_interval.first; i < out_interval.second; i++){
               uint64_t u = out_e[i];
               histogram[labels0[u]]++;
           }

           // cfr. Spec v0.9 pp 14 "If the graph is directed and a neighbor is reachable via both an incoming and
           // outgoing edge, its label will be counted twice"
           if(m_is_directed){
               auto in_interval = get_in_interval(v);
               for(uint64_t i = in_interval.first; i < in_interval.second; i++){
                   uint64_t u = in_e[i];
                   histogram[labels0[u]]++;
               }
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

void CSR::cdlp(uint64_t max_iterations, const char* dump2file) {
    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();

    // Run the CDLP algorithm
    unique_ptr<uint64_t[]> labels = do_cdlp(max_iterations, timeout);
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);  }

    // Translate the vertex IDs
    auto translation = translate(labels.get(), m_num_vertices);
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // Store the results in the given file
    if(dump2file != nullptr){
        save_results(translation, dump2file);
    }
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

unique_ptr<double[]> CSR::do_lcc(utility::TimeoutService& timer) const {
    if(m_is_directed){
        return do_lcc_directed(timer);
    } else {
        return do_lcc_undirected(timer);
    }
}

// loosely based on the impl~ made for GraphOne
unique_ptr<double[]> CSR::do_lcc_directed(utility::TimeoutService& timer) const {
    assert(m_is_directed && "Implementation for directed graphs");

    unique_ptr<double[]> ptr_lcc { new double[m_num_vertices] };
    double* lcc = ptr_lcc.get();
    uint64_t* __restrict out_e = m_out_e;
    uint64_t* __restrict in_e = m_in_e;

    #pragma omp parallel
    {
        std::vector<uint64_t> edges;

        #pragma omp for schedule(dynamic, 64)
        for(uint64_t v = 0; v < m_num_vertices; v++){
            COUT_DEBUG_LCC("> Node " << v);
            if(timer.is_timeout()) continue; // exhausted the budget of available time
            lcc[v] = 0.0;
            uint64_t num_triangles = 0; // number of triangles found so far for the node v

            // Cfr. Spec v.0.9.0 pp. 15: "If the number of neighbors of a vertex is less than two, its coefficient is defined as zero"
            uint64_t v_degree_out = get_out_degree(v);
            uint64_t v_degree_in = get_in_degree(v);
            uint64_t v_degree_ub = v_degree_in + v_degree_out; // upper bound for directed graphs, exact degree for those undirected
            if(v_degree_ub < 2) continue;

            // Build the list of neighbours of v
            unordered_set<uint64_t> neighbours;
            edges.clear();
            edges.reserve(v_degree_ub);

            // Outgoing edges
            auto out_interval = get_out_interval(v);
            for(uint64_t i = out_interval.first; i < out_interval.second; i++){
                uint64_t u = out_e[i];
                edges.push_back(u);
                neighbours.insert(u);
            }

            // Incoming edges (only directed graphs)
            auto in_interval = get_in_interval(v);
            for(uint64_t i = in_interval.first; i < in_interval.second; i++){
                uint64_t u = in_e[i];
                auto result = neighbours.insert(u);
                if(result.second){ // the element was actually inserted
                    edges.push_back(u);
                }
            }
            const uint64_t v_degree = edges.size();

            // Now we know is the actual degree of v, perform the proper check for directed graphs
            if(v_degree  < 2) continue;

            // again, visit all neighbours of v
            // for directed graphs, edges1 contains the intersection of both the incoming and the outgoing edges
            for(uint64_t i = 0, end = v_degree; i < end; i++){
                uint64_t u = edges[i];
                COUT_DEBUG_LCC("[" << i << "/" << edges.size() << "] neighbour: " << u);
                assert(neighbours.count(u) == 1 && "The set `neighbours' should contain all neighbours of v");

                // For the Graphalytics spec v 0.9.0, only consider the outgoing edges for the neighbours u
                auto u_out_interval = get_out_interval(u);

                for(uint64_t j = u_out_interval.first; j < u_out_interval.second; j++){
                    uint64_t w = out_e[j];
                    COUT_DEBUG_LCC("---> [" << j << "/" << /* degree */ (u_out_interval.second - u_out_interval.first) << "] neighbour: " << w);
                    // check whether it's also a neighbour of v
                    if(neighbours.count(w) == 1){
                        COUT_DEBUG_LCC("Triangle found " << v << " - " << u << " - " << w);
                        num_triangles++;
                    }
                }
            }

            // register the final score
            uint64_t max_num_edges = v_degree * (v_degree -1);
            lcc[v] = static_cast<double>(num_triangles) / max_num_edges;
            COUT_DEBUG_LCC("Score computed: " << (num_triangles) << "/" << max_num_edges << " = " << lcc[v]);
        }
    }

    return ptr_lcc;
}

unique_ptr<double[]> CSR::do_lcc_undirected(utility::TimeoutService& timer) const {
    assert(!m_is_directed && "Implementation for undirected graphs");

    unique_ptr<double[]> ptr_lcc { new double[m_num_vertices] };
    double* lcc = ptr_lcc.get();
    uint64_t* __restrict out_e = m_out_e;

    #pragma omp parallel for schedule(dynamic, 64)
    for(uint64_t v = 0; v < m_num_vertices; v++){
        COUT_DEBUG_LCC("> Node " << v);
        if(timer.is_timeout()) continue; // exhausted the budget of available time
        lcc[v] = 0.0;
        uint64_t num_triangles = 0; // number of triangles found so far for the node v

        // Cfr. Spec v.0.9.0 pp. 15: "If the number of neighbors of a vertex is less than two, its coefficient is defined as zero"
        uint64_t v_degree_out = get_out_degree(v);
        if(v_degree_out < 2) continue;

        // Build the list of neighbours of v
        unordered_set<uint64_t> neighbours;

        // Outgoing edges
        auto out_interval = get_out_interval(v);
        for(uint64_t i = out_interval.first; i < out_interval.second; i++){
            uint64_t u = out_e[i];
            neighbours.insert(u);
        }

        // again, visit all neighbours of v
        for(uint64_t i = out_interval.first; i < out_interval.second; i++){
            uint64_t u = out_e[i];
            COUT_DEBUG_LCC("[" << (i - out_interval.first) << "/" << v_degree_out << "] neighbour: " << u);
            assert(neighbours.count(u) == 1 && "The set `neighbours' should contain all neighbours of v");

            // For the Graphalytics spec v 0.9.0, only consider the outgoing edges for the neighbours u
            auto u_out_interval = get_out_interval(u);

            for(uint64_t j = u_out_interval.first; j < u_out_interval.second; j++){
                uint64_t w = out_e[j];
                COUT_DEBUG_LCC("---> [" << j << "/" << /* degree */ (u_out_interval.second - u_out_interval.first) << "] neighbour: " << w);
                // check whether it's also a neighbour of v
                if(neighbours.count(w) == 1){
                    COUT_DEBUG_LCC("Triangle found " << v << " - " << u << " - " << w);
                    num_triangles++;
                }
            }
        }

        // register the final score
        uint64_t max_num_edges = v_degree_out * (v_degree_out -1);
        lcc[v] = static_cast<double>(num_triangles) / max_num_edges;
        COUT_DEBUG_LCC("Score computed: " << (num_triangles) << "/" << max_num_edges << " = " << lcc[v]);
    }
    return ptr_lcc;
}

void CSR::lcc(const char* dump2file) {
    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();

    // Run the LCC algorithm
    unique_ptr<double[]> scores = do_lcc(timeout);
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);  }

    auto translation = translate(scores.get(), m_num_vertices);
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // Store the results in the given file
    if(dump2file != nullptr){
        save_results(translation, dump2file);
    }
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

gapbs::pvector<WeightT> CSR::do_sssp(uint64_t source, double delta, utility::TimeoutService& timer) const {
    // Init
    gapbs::pvector<WeightT> dist(num_vertices(), numeric_limits<WeightT>::infinity());
    dist[source] = 0;
    gapbs::pvector<NodeID> frontier(num_edges());
    // two element arrays for double buffering curr=iter&1, next=(iter+1)&1
    size_t shared_indexes[2] = {0, kMaxBin};
    size_t frontier_tails[2] = {1, 0};
    frontier[0] = source;
    uint64_t* __restrict out_e = m_out_e;
    double* __restrict out_w = m_out_w;

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
                    const auto u_interval = get_out_interval(u);
                    for(uint64_t i = u_interval.first; i < u_interval.second; i++){
                        uint64_t v = out_e[i];
                        double w = out_w[i];

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

void CSR::sssp(uint64_t source_vertex_id, const char* dump2file) {
    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();

    // Run the SSSP algorithm
    double delta = 2.0; // same value used in the GAPBS, at least for most graphs
    auto distances = do_sssp(m_ext2log.at(source_vertex_id), delta, timeout);
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);  }

    // Translate the logical IDs into the external IDs
    auto translation = translate(distances.data(), m_num_vertices);
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // store the results in the given file
    if(dump2file != nullptr)
        save_results(translation, dump2file);
}


/*****************************************************************************
 *                                                                           *
 *  LCC, sort-merge implementation                                           *
 *                                                                           *
 *****************************************************************************/
/**
 * Algorithm parameters
 */
static const uint64_t LCC_NUM_WORKERS = thread::hardware_concurrency(); // number of workers / logical threads to use
//static const uint64_t LCC_NUM_WORKERS = 1; // number of workers / logical threads to use
static constexpr uint64_t LCC_TASK_SIZE = 1ull << 10; // number of vertices processed in each task

#undef COUT_CLASS_NAME
#define COUT_CLASS_NAME "CSR_LCC"

CSR_LCC::CSR_LCC(bool is_directed, bool numa_interleaved) : CSR(is_directed, numa_interleaved) {
    if(is_directed) { ERROR("[CSR_LCC] Directed graphs are not supported"); }
}

void CSR_LCC::lcc(const char* dump2file){
    // Init
    unique_ptr<double[]> scores;
    utility::TimeoutService timeout { m_timeout };
    common::Timer timer; timer.start();

    { // restrict the scope to allow the dtor to clean up
        Master algorithm ( this, &timeout );
        scores.reset( algorithm.execute() );
    }
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);  }

    // Translate the vertex IDs
    auto translation = translate(scores.get(), m_num_vertices);
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // store the results in the given file
    if(dump2file != nullptr){
        save_results(translation, dump2file);
    }
}

#undef COUT_CLASS_NAME
#define COUT_CLASS_NAME "CSR_LCC::Master"

CSR_LCC::Master::Master(const CSR_LCC* csr, utility::TimeoutService* timeout) :
        m_csr(csr), m_scores(nullptr), m_num_triangles(nullptr), m_next(0), m_timeout(timeout){
}

CSR_LCC::Master::~Master(){
    delete[] m_scores; m_scores = nullptr;
    delete[] m_num_triangles; m_num_triangles = nullptr;
}

void CSR_LCC::Master::initialise(){
    assert(m_num_triangles == nullptr && "Already initialised");

    m_scores = new double[m_csr->num_vertices()](); // init to 0
    m_num_triangles = new atomic<uint64_t>[m_csr->num_vertices()](); // init to 0;
}

void CSR_LCC::Master::compute_scores(){
    //common::Timer timer; timer.start();

    for(uint64_t i = 0, N = m_csr->num_vertices(); i < N; i++){
        uint64_t num_triangles = m_num_triangles[i];
        if(num_triangles > 0){
            uint64_t degree = m_csr->get_out_degree(i);
            uint64_t max_num_edges = degree * (degree -1);
            double score = static_cast<double>(num_triangles) / max_num_edges;
            COUT_DEBUG("vertex: " << i << ", external id: " << m_csr->m_log2ext[i] << ", num triangles: " << num_triangles << ", degree: " << degree << ", score: " << score);
            m_scores[i] = score;
        } // else m_scores[i] = 0 (default value)
    }

    //timer.stop();
    //COUT_DEBUG_FORCE("compute_scores executed in: " << timer);
}

double* CSR_LCC::Master::execute() {
    common::Timer timer; timer.start();

    // init the state and the side information for each vertex
    initialise();
    //COUT_DEBUG_FORCE("Initialisation time: " << timer);

    // start the workers
    assert(LCC_NUM_WORKERS >= 1 && "At least one worker should be set");
    vector<Worker*> workers;
    workers.reserve(LCC_NUM_WORKERS);
    for(uint64_t worker_id = 0; worker_id < LCC_NUM_WORKERS; worker_id++ ){
        workers.push_back(new Worker(m_csr, this));
    }

    // wait for the workers to terminate ...
    for(uint64_t worker_id = 0; worker_id < workers.size(); worker_id++ ){
        workers[worker_id]->join();
        delete workers[worker_id]; workers[worker_id] = nullptr;
    }
    if(m_timeout->is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    compute_scores();

    double* scores = m_scores;
    m_scores = nullptr;
    return scores;
}

bool CSR_LCC::Master::next_task(uint64_t* output_vtx_start /* inclusive */, uint64_t* output_vtx_end /* exclusive */) {
    uint64_t logical_start = m_next.fetch_add(LCC_TASK_SIZE); /* return the previous value of m_next */
    uint64_t num_vertices = m_csr->num_vertices();
    if(logical_start >= num_vertices || m_timeout->is_timeout()){
        return false;
    } else {
        uint64_t logical_end = std::min(logical_start + LCC_TASK_SIZE, num_vertices);

        *output_vtx_start = logical_start;
        *output_vtx_end = logical_end;

        return true;
    }
}

atomic<uint64_t>& CSR_LCC::Master::num_triangles(uint64_t vertex_id) {
    return m_num_triangles[vertex_id];
}

#undef COUT_CLASS_NAME
#define COUT_CLASS_NAME "CSR_LCC::Worker"

CSR_LCC::Worker::Worker(const CSR_LCC* csr, Master* master) : m_csr(csr), m_master(master) {
    m_handle = thread { &Worker::execute, this };
}

CSR_LCC::Worker::~Worker(){ }

void CSR_LCC::Worker::execute() {
    COUT_DEBUG("Worker started");

    uint64_t v_start, v_end;
    while(m_master->next_task(&v_start, &v_end)){
        for(uint64_t v = v_start; v < v_end; v++){
            process_vertex(v);
        }
    }

    COUT_DEBUG("Worker terminated");
}

void CSR_LCC::Worker::join(){
    m_handle.join();
}

void CSR_LCC::Worker::process_vertex(uint64_t n1) {
    COUT_DEBUG("vertex: " << n1);
    uint64_t num_triangles = 0; // current number of triangles found for `n1'
    m_neighbours.clear();
    uint64_t* __restrict out_e = m_csr->m_out_e;

    auto n1_interval = m_csr->get_out_interval(n1);
    for(uint64_t i = n1_interval.first; i < n1_interval.second; i++){
        uint64_t n2 = out_e[i];
        if(n2 > n1) break; // we're done with n1

        m_neighbours.push_back(n2);
        uint64_t marker = 0; // current position in the neighbours vector, to merge shared neighbours

        auto n2_interval = m_csr->get_out_interval(n2);
        for(uint64_t j = n2_interval.first; j < n2_interval.second; j++){
            uint64_t n3 = out_e[j];
            if(n3 > n2) break; // we're done with n2
            assert(n1 > n2 && n2 > n3); // we're looking for triangles of the kind c - b - a, with c > b && b > a

            COUT_DEBUG("  candidate: " << n1 << " - " << n2 << " - " << n3 << ", marker[" << marker << "] = " << m_neighbours[marker]);

            if(n3 > m_neighbours[marker]){ // merge with m_neighbours
                do {
                    marker ++;
                } while(marker < m_neighbours.size() && n3 > m_neighbours[marker]);
                if(marker >= m_neighbours.size()) break; // there is nothing left to merge
            }

            if(n3 == m_neighbours[marker]){ // match !
                COUT_DEBUG("    match: " << n1 << " - " << n2 << " - " << n3 << " and " << n1 << " - " << n3 << " - " << n2);

                num_triangles += 2; // we've discovered both n1 - n2 - n3 and n1 - n3 - n2; with n1 > n2 > n3

                // increase the contribution for n2
                m_master->num_triangles(n2) += 2;

                // increase the contribution for n3
                m_master->num_triangles(n3) += 2;

                marker++;
                if(marker >= m_neighbours.size()) break; // there is nothing left to merge
            }
        }
    }

    if(num_triangles != 0){
        m_master->num_triangles(n1) += num_triangles;
    }
}

} // namespace


