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
#include "teseo_driver.hpp"

#include <cassert>
#include <fstream>
#include <iomanip>
#include <mutex>
#include <omp.h>
#include <unordered_map>
#include <unordered_set>

#include "common/timer.hpp"
#include "third-party/gapbs/gapbs.hpp"
#include "third-party/libcuckoo/cuckoohash_map.hh"
#include "utility/timeout_service.hpp"
#include "teseo_openmp.hpp"
#include "teseo/context/global_context.hpp"
#include "teseo/memstore/memstore.hpp"
#include "teseo.hpp"

using namespace common;
using namespace gapbs;
using namespace libcuckoo;
using namespace std;
using namespace teseo;
using namespace gfe::library::teseo_driver_internal;

#define TESEO reinterpret_cast<Teseo*>(m_pImpl)

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG
namespace gfe { extern mutex _log_mutex [[maybe_unused]]; }
#define COUT_DEBUG_FORCE(msg) { std::scoped_lock<std::mutex> lock{::gfe::_log_mutex}; std::cout << "[TeseoDriver::" << __FUNCTION__ << "] " << msg << std::endl; }
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
TeseoDriver::TeseoDriver(bool is_directed, bool read_only) : m_pImpl(new Teseo()), m_is_directed(is_directed), m_read_only(read_only), m_thread_affinity(true) {
    if(is_directed == true){ throw std::invalid_argument("Only undirected graphs are currently supported by the front-end"); }
}

TeseoDriver::~TeseoDriver(){
    delete TESEO; m_pImpl = nullptr;
}

void TeseoDriver::on_thread_init(int thread_id){
    TESEO->register_thread();
}

void TeseoDriver::on_thread_destroy(int thread_id){
    TESEO->unregister_thread();
}

void* TeseoDriver::handle_impl(){
    return m_pImpl;
}

/*****************************************************************************
 *                                                                           *
 *  Updates & point look ups                                                 *
 *                                                                           *
 *****************************************************************************/

uint64_t TeseoDriver::num_edges() const {
    auto tx = TESEO->start_transaction();
    return tx.num_edges();
}

uint64_t TeseoDriver::num_vertices() const {
    auto tx = TESEO->start_transaction();
    return tx.num_vertices();
}

bool TeseoDriver::has_vertex(uint64_t vertex_id) const {
    auto tx = TESEO->start_transaction();
    return tx.has_vertex(vertex_id);
}

double TeseoDriver::get_weight(uint64_t source, uint64_t destination) const {
    auto tx = TESEO->start_transaction();
    try {
        return tx.get_weight(source, destination);
    } catch (LogicalError& e) {
        return numeric_limits<double>::signaling_NaN();
    }
}

bool TeseoDriver::is_directed() const {
    return m_is_directed;
}

void TeseoDriver::set_timeout(uint64_t seconds){
    m_timeout = chrono::seconds{ seconds };
}

void TeseoDriver::set_thread_affinity(bool value){
    m_thread_affinity = value;
}

bool TeseoDriver::has_thread_affinity() const {
    return m_thread_affinity;
}

bool TeseoDriver::has_read_only_transactions() const {
    return m_read_only;
}

bool TeseoDriver::add_vertex(uint64_t vertex_id){
    auto tx = TESEO->start_transaction();
    try {
        tx.insert_vertex(vertex_id);
        tx.commit();
        return true;
    } catch( LogicalError& e ){
        return false;
    } catch( TransactionConflict& e) {
        return false;
    }
}

bool TeseoDriver::remove_vertex(uint64_t vertex_id){
    COUT_DEBUG("remove vertex: " << vertex_id);
    auto tx = TESEO->start_transaction();
    try {
        tx.remove_vertex(vertex_id);
        tx.commit();
        return true;
    } catch( LogicalError& e ){
        return false;
    } catch( TransactionConflict& e) {
        return false;
    }
}

bool TeseoDriver::add_edge(gfe::graph::WeightedEdge e) {
    if(e.source() == e.destination())
        throw std::invalid_argument("Teseo does not support self-loop edges (e.g. same source, same destination)");

    auto tx = TESEO->start_transaction();
    try {
        tx.insert_edge(e.source(), e.destination(), e.weight());
        tx.commit();
        return true;
    } catch( LogicalError& e ){
        return false;
    } catch( TransactionConflict& e) {
        return false;
    }
}

bool TeseoDriver::add_edge_v2(gfe::graph::WeightedEdge e){
    if(e.source() == e.destination())
        throw std::invalid_argument("Teseo does not support self-loop edges (e.g. same source, same destination)");

    static cuckoohash_map<uint64_t, bool> vertices;

    if(vertices.insert(e.source(), true)){ add_vertex(e.source()); }
    if(vertices.insert(e.destination(), true)){ add_vertex(e.destination()); }
    while( ! add_edge(e) ) { /* nop */ }

    return true;
}

bool TeseoDriver::remove_edge(gfe::graph::Edge e) {
    auto tx = TESEO->start_transaction();
    try {
        tx.remove_edge(e.source(), e.destination());
        tx.commit();
        return true;
    } catch( LogicalError& e ){
        return false;
    } catch( TransactionConflict& e) {
        return false;
    }
}

/*****************************************************************************
 *                                                                           *
 *  Dump                                                                     *
 *                                                                           *
 *****************************************************************************/

void TeseoDriver::dump_ostream(std::ostream& out) const {
    auto memstore = teseo::context::global_context()->memstore();
    memstore->dump();
}


/*****************************************************************************
 *                                                                           *
 *  Helpers                                                                  *
 *                                                                           *
 *****************************************************************************/
// Translate the logical into real vertices IDs. Materialization step at the end of a Graphalytics algorithm
template <typename T>
static vector<pair<uint64_t, T>> translate(OpenMP& openmp, T* values, uint64_t N){
    vector<pair<uint64_t , T>> external_ids(N);

    #pragma omp parallel for firstprivate(openmp)
    for (uint64_t v = 0; v < N; v++) {
        external_ids[v] = make_pair(openmp.transaction().vertex_id(v), values[v]);
    }

    return external_ids;
}

template <typename T, bool negative_scores = true>
void TeseoDriver::save_results(std::vector<std::pair<uint64_t, T>>& result, const char* dump2file) {
    if(dump2file == nullptr) return; // nop

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
namespace { // anonymous

/*
GAP Benchmark Suite
Kernel: Breadth-First Search (BFS)
Author: Scott Beamer
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

static
int64_t BUStep(OpenMP& openmp, pvector<int64_t>& distances, int64_t distance, Bitmap &front, Bitmap &next) {

    const int64_t N = openmp.transaction().num_vertices();
    int64_t awake_count = 0;
    next.reset();
    #pragma omp parallel for reduction(+ : awake_count) schedule(dynamic, 1024) firstprivate(openmp)
    for (int64_t u=0; u < N; u++) {
        if (distances[u] < 0){ // the node has not been visited yet

            bool done = false;
            openmp.iterator().edges(u, true, [u, &done, &awake_count, &distances, distance, &front, &next](uint64_t n){

                if (front.get_bit(n)) {
                    distances[u] = distance; // on each BUStep, all nodes will have the same distance
                    awake_count++;
                    next.set_bit(u);
                    done = true;
                }

                return !done;
            });
        }
    }
    return awake_count;
}

static
int64_t TDStep(OpenMP& openmp, pvector<int64_t>& distances, int64_t distance, SlidingQueue<int64_t>& queue) {
    int64_t scout_count = 0;

    #pragma omp parallel firstprivate(openmp)
    {
        QueueBuffer<int64_t> lqueue(queue);
        #pragma omp for reduction(+ : scout_count)
        for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
            int64_t u = *q_iter;

            openmp.iterator().edges(u, true, [&distances, distance, &lqueue, &scout_count](uint64_t destination){
                int64_t curr_val = distances[destination];

                if (curr_val < 0 && compare_and_swap(distances[destination], curr_val, distance)) {
                    lqueue.push_back(destination);
                    scout_count += -curr_val;
                }
            });
        }
        lqueue.flush();
    }

    return scout_count;
}

static
void QueueToBitmap(const SlidingQueue<int64_t> &queue, Bitmap &bm) {
    #pragma omp parallel for
    for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
        int64_t u = *q_iter;
        bm.set_bit_atomic(u);
    }
}

static
void BitmapToQueue(OpenMP& openmp, const Bitmap &bm, SlidingQueue<int64_t> &queue) {
    const int64_t N = openmp.transaction().num_vertices();

    #pragma omp parallel
    {
        QueueBuffer<int64_t> lqueue(queue);
        #pragma omp for
        for (int64_t n=0; n < N; n++)
            if (bm.get_bit(n))
                lqueue.push_back(n);
        lqueue.flush();
    }
    queue.slide_window();
}

static
pvector<int64_t> InitDistances(OpenMP& openmp){
    const int64_t N = openmp.transaction().num_vertices();
    pvector<int64_t> distances(N);
    #pragma omp parallel for firstprivate(openmp)
    for (int64_t n = 0; n < N; n++){
        int64_t out_degree = openmp.transaction().degree(n, /* logical ? */ true);
        distances[n] = out_degree != 0 ? - out_degree : -1;
    }
    return distances;
}

} // anon namespace

static
pvector<int64_t> teseo_bfs(OpenMP& openmp, int64_t source, utility::TimeoutService& timer, int alpha = 15, int beta = 18) {
    // The implementation from GAP BS reports the parent (which indeed it should make more sense), while the one required by
    // Graphalytics only returns the distance

    pvector<int64_t> distances = InitDistances(openmp);
    distances[source] = 0;

    SlidingQueue<int64_t> queue(openmp.transaction().num_vertices());
    queue.push_back(source);
    queue.slide_window();
    Bitmap curr(openmp.transaction().num_vertices());
    curr.reset();
    Bitmap front(openmp.transaction().num_vertices());
    front.reset();
    int64_t edges_to_check = openmp.transaction().num_edges();
    int64_t scout_count = openmp.transaction().degree(source, true);
    int64_t distance = 1; // current distance
    while (!timer.is_timeout() && !queue.empty()) {

        if (scout_count > edges_to_check / alpha) {
            int64_t awake_count, old_awake_count;
            QueueToBitmap(queue, front);
            awake_count = queue.size();
            queue.slide_window();
            do {
                old_awake_count = awake_count;
                awake_count = BUStep(openmp, distances, distance, front, curr);
                front.swap(curr);
                distance++;
            } while ((awake_count >= old_awake_count) || (awake_count > static_cast<int64_t>(openmp.transaction().num_vertices()) / beta));
            BitmapToQueue(openmp, front, queue);
            scout_count = 1;
        } else {
            edges_to_check -= scout_count;
            scout_count = TDStep(openmp, distances, distance, queue);
            queue.slide_window();
            distance++;
        }
    }

    return distances;
}

void TeseoDriver::bfs(uint64_t source_vertex_id, const char* dump2file){
    OpenMP openmp ( this );

    utility::TimeoutService tcheck { m_timeout };
    common::Timer timer; timer.start();

    // execute the BFS algorithm
    auto result = teseo_bfs(openmp, openmp.transaction().logical_id(source_vertex_id), tcheck);
    if(tcheck.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // translate the logical IDs into the external IDs
    const int64_t N = openmp.transaction().num_vertices();
    vector<pair<uint64_t , int64_t>> external_ids(N);
    #pragma omp parallel for firstprivate(openmp)
    for (int64_t v = 0; v < N; v++) {
        external_ids[v] = make_pair(openmp.transaction().vertex_id(v), result[v]);
    }

    // store the results in the given file
    if(dump2file != nullptr)
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
it is not necesarily the fastest way to implement it. It does perform the
updates in the pull direction to remove the need for atomics.
*/

static
unique_ptr<double[]> teseo_pagerank(OpenMP& openmp, uint64_t num_iterations, double damping_factor, utility::TimeoutService& timer) {
    // init
    const uint64_t num_vertices = openmp.transaction().num_vertices();
    COUT_DEBUG("num vertices: " << num_vertices);

    const double init_score = 1.0 / num_vertices;
    const double base_score = (1.0 - damping_factor) / num_vertices;

    unique_ptr<double[]> ptr_scores{ new double[num_vertices]() }; // avoid memory leaks
    double* scores = ptr_scores.get();
    #pragma omp parallel for
    for(uint64_t v = 0; v < num_vertices; v++){
        scores[v] = init_score;
    }
    gapbs::pvector<double> outgoing_contrib(num_vertices, 0.0);

    // pagerank iterations
    for(uint64_t iteration = 0; iteration < num_iterations && !timer.is_timeout(); iteration++){
        double dangling_sum = 0.0;

        // for each node, precompute its contribution to all of its outgoing neighbours and, if it's a sink,
        // add its rank to the `dangling sum' (to be added to all nodes).
        #pragma omp parallel for reduction(+:dangling_sum) firstprivate(openmp)
        for(uint64_t v = 0; v < num_vertices; v++){
            uint64_t out_degree = openmp.transaction().degree(v, /* logical */ true);
            if(out_degree == 0){ // this is a sink
                dangling_sum += scores[v];
            } else {
                outgoing_contrib[v] = scores[v] / out_degree;
            }
        }

        dangling_sum /= num_vertices;

        // compute the new score for each node in the graph
        #pragma omp parallel for schedule(dynamic, 64) firstprivate(openmp)
        for(uint64_t v = 0; v < num_vertices; v++){

            double incoming_total = 0;
            openmp.iterator().edges(v, /* logical ? */ true, [&incoming_total, &outgoing_contrib](uint64_t destination){
               incoming_total += outgoing_contrib[destination];
            });

            // update the score
            scores[v] = base_score + damping_factor * (incoming_total + dangling_sum);
        }
    }

    return ptr_scores;
}

void TeseoDriver::pagerank(uint64_t num_iterations, double damping_factor, const char* dump2file) {
    OpenMP openmp( this ); // OpenMP machinery

    // Init
    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();

    // Run the PageRank algorithm
    unique_ptr<double[]> ptr_rank = teseo_pagerank(openmp, num_iterations, damping_factor, timeout);
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);  }

    // translate the vertex IDs
    const uint64_t N = openmp.transaction().num_vertices();
    auto external_ids = translate(openmp, ptr_rank.get(), N);
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // store the results in the given file
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
static // do_wcc
unique_ptr<uint64_t[]> teseo_wcc(OpenMP& openmp, utility::TimeoutService& timer) {
    // init
    const uint64_t N = openmp.transaction().num_vertices();
    unique_ptr<uint64_t[]> ptr_components { new uint64_t[N] };
    uint64_t* comp = ptr_components.get();

    #pragma omp parallel for
    for (uint64_t n = 0; n < N; n++){
        comp[n] = n;
    }

    bool change = true;
    while (change && !timer.is_timeout()) {
        change = false;

        #pragma omp parallel for firstprivate(openmp)
        for (uint64_t u = 0; u < N; u++){

            openmp.iterator().edges(u, true, [comp, u, &change](uint64_t v){
                uint64_t comp_u = comp[u];
                uint64_t comp_v = comp[v];
                if (comp_u == comp_v) return true;

                // Hooking condition so lower component ID wins independent of direction
                uint64_t high_comp = std::max(comp_u, comp_v);
                uint64_t low_comp = std::min(comp_u, comp_v);
                if (high_comp == comp[high_comp]) {
                    change = true;
                    comp[high_comp] = low_comp;
                }

                return true;
            });

        }

        #pragma omp parallel for
        for (uint64_t n = 0; n < N; n++){
            while (comp[n] != comp[comp[n]]) {
                comp[n] = comp[comp[n]];
            }
        }
    }

    return ptr_components;
}

void TeseoDriver::wcc(const char* dump2file) {
    OpenMP openmp ( this ); // OpenMP machinery

    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();

    // run wcc
    unique_ptr<uint64_t[]> ptr_components = teseo_wcc(openmp, timeout);

    // translate the vertex IDs
    const uint64_t N = openmp.transaction().num_vertices();
    auto external_ids = translate<uint64_t>(openmp, ptr_components.get(), N);
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
static unique_ptr<uint64_t[]> teseo_cdlp(OpenMP& openmp, uint64_t max_iterations, utility::TimeoutService& timer){
    const uint64_t num_vertices = openmp.transaction().num_vertices();
    unique_ptr<uint64_t[]> ptr_labels0 { new uint64_t[num_vertices] };
    unique_ptr<uint64_t[]> ptr_labels1 { new uint64_t[num_vertices] };
    uint64_t* labels0 = ptr_labels0.get(); // current labels
    uint64_t* labels1 = ptr_labels1.get(); // labels for the next iteration

    // init
    #pragma omp parallel for firstprivate(openmp)
    for(uint64_t v = 0; v < num_vertices; v++){
        labels0[v] = openmp.transaction().vertex_id(v);
    }

    // algorithm pass
    bool change = true;
    uint64_t current_iteration = 0;
    while(current_iteration < max_iterations && change && !timer.is_timeout()){
        change = false; // reset the flag

        #pragma omp parallel for shared(change) firstprivate(openmp)
        for(uint64_t v = 0; v < num_vertices; v++){
            unordered_map<uint64_t, uint64_t> histogram;

            // compute the histogram from both the outgoing & incoming edges. The aim is to find the number of each label
            // shared among the neighbours of node_id
            openmp.iterator().edges(v, true, [&histogram, labels0](uint64_t u){
                histogram[labels0[u]]++;
            });

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

void TeseoDriver::cdlp(uint64_t max_iterations, const char* dump2file) {
    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();

    OpenMP openmp ( this ); // OpenMP thread-state, init the transaction

    // Run the CDLP algorithm
    unique_ptr<uint64_t[]> labels = teseo_cdlp(openmp, max_iterations, timeout);
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);  }

    // translate the vertex IDs
    const uint64_t N = openmp.transaction().num_vertices();
    auto external_ids = translate<uint64_t>(openmp, labels.get(), N);
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // store the results in the given file
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

// loosely based on the impl~ made for Stinger
static unique_ptr<double[]> teseo_lcc(OpenMP& openmp, utility::TimeoutService& timer){
    const uint64_t num_vertices = openmp.transaction().num_vertices();

    unique_ptr<double[]> ptr_lcc { new double[num_vertices] };
    double* lcc = ptr_lcc.get();

    #pragma omp parallel for firstprivate(openmp)
    for(uint64_t v = 0; v < num_vertices; v++){
        COUT_DEBUG_LCC("> Node " << v);
        if(timer.is_timeout()) continue; // exhausted the budget of available time

        lcc[v] = 0.0;
        uint64_t num_triangles = 0; // number of triangles found so far for the node v

        // Cfr. Spec v.0.9.0 pp. 15: "If the number of neighbors of a vertex is less than two, its coefficient is defined as zero"
        uint64_t degree = openmp.transaction().degree(v, true);
        if(degree < 2) continue;

        // Build the list of neighbours of v
        unordered_set<uint64_t> neighbours;

        openmp.iterator().edges(v, true, [&neighbours](uint64_t destination){
            neighbours.insert(destination);
        });

        // again, visit all neighbours of v
        openmp.iterator().edges(v, true, [&neighbours, &num_triangles, &openmp](uint64_t u){
            assert(neighbours.count(u) == 1 && "The set `neighbours' should contain all neighbours of v");

            openmp.iterator().edges(u, true, [&neighbours, &num_triangles](uint64_t w){
                // check whether it's also a neighbour of v
                if(neighbours.count(w) == 1){
                    //COUT_DEBUG_LCC("Triangle found " << v << " - " << u << " - " << w);
                    //COUT_DEBUG_LCC("  -> " << transaction.vertex_id(v) << " - " << transaction.vertex_id(u) << " - " << transaction.vertex_id(w));
                    num_triangles++;
                }
            });
        });

        // register the final score
        uint64_t max_num_edges = degree * (degree -1);
        lcc[v] = static_cast<double>(num_triangles) / max_num_edges;
        COUT_DEBUG_LCC("Score computed: " << (num_triangles) << "/" << max_num_edges << " = " << lcc[v]);
    }

    return ptr_lcc;
}

void TeseoDriver::lcc(const char* dump2file) {
    OpenMP openmp { this }; // OpenMP machinery

    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();

    // Run the LCC algorithm
    unique_ptr<double[]> scores = teseo_lcc(openmp, timeout);
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);  }

    // translate the vertex IDs
    const uint64_t N = openmp.transaction().num_vertices();
    auto external_ids = translate<double>(openmp, scores.get(), N);
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // store the results in the given file
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

static gapbs::pvector<WeightT> teseo_sssp(OpenMP& openmp, uint64_t source, double delta, utility::TimeoutService& timer){
    const uint64_t num_vertices = openmp.transaction().num_vertices();
    const uint64_t num_edges = openmp.transaction().num_edges();

    // Init
    gapbs::pvector<WeightT> dist(num_vertices, numeric_limits<WeightT>::infinity());
    dist[source] = 0;
    gapbs::pvector<NodeID> frontier(num_edges);
    // two element arrays for double buffering curr=iter&1, next=(iter+1)&1
    size_t shared_indexes[2] = {0, kMaxBin};
    size_t frontier_tails[2] = {1, 0};
    frontier[0] = source;

    #pragma omp parallel firstprivate(openmp)
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
                    openmp.iterator().edges(u, /* logical ? */ true, [u, delta, &dist, &local_bins](uint64_t v, double w){
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


                        return true;
                    });
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

void TeseoDriver::sssp(uint64_t source_vertex_id, const char* dump2file) {
    OpenMP openmp { this }; // OpenMP machinery, init the transaction

    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();

    // Run the SSSP algorithm
    double delta = 2.0; // same value used in the GAPBS, at least for most graphs
    auto distances = teseo_sssp(openmp, openmp.transaction().logical_id(source_vertex_id), delta, timeout);
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);  }

    // translate the vertex IDs
    const uint64_t N = openmp.transaction().num_vertices();
    vector<pair<uint64_t, double>> external_ids(N);
    #pragma omp parallel for firstprivate(openmp)
    for (uint64_t v = 0; v < N; v++) {
        external_ids[v] = make_pair(openmp.transaction().vertex_id(v), distances[v]);
    }
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // store the results in the given file
    if(dump2file != nullptr)
        save_results(external_ids, dump2file);
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

namespace {

/**
 * Forward declarations
 */
class LCC_Master;
class LCC_Worker;

class LCC_Master {
    Teseo * const m_teseo; // instance to the database
    Transaction m_transaction; // the transaction providing isolation for the computation
    double* m_scores; // the final scores of the LCC algorithm
    std::atomic<uint64_t>* m_num_triangles; // number of triangles counted so far for the given vertex, array of num_vertices
    std::atomic<uint64_t> m_next; // counter to select the next task among the workers
    utility::TimeoutService* m_timeout; // timer to check whether we are not spending more time than what allocated (1 hour typically)

    // Reserve the space in the hash maps m_score and m_state so that they can be operated concurrently by each thread/worker
    void initialise();

    // Compute the final scores
    void compute_scores();

public:
    // Constructor
    LCC_Master(Teseo* teseo, Transaction transaction, utility::TimeoutService* timeout);

    // Destructor
    ~LCC_Master();

    // Execute the algorithm
    double* execute();

    // Select the next window to process, in the form [vertex_start, vertex_end);
    // Return true if a window/task has been fetched, false if there are no more tasks to process
    bool next_task(uint64_t* output_vtx_start /* inclusive */, uint64_t* output_vtx_end /* exclusive */);

    // Retrieve the instance to the database
    Teseo* teseo();

    // Retrieve the number of triangles associated to the given vertex
    std::atomic<uint64_t>& num_triangles(uint64_t vertex_id);
};

class LCC_Worker {
    LCC_Master* m_master; // handle to the master instance
    Transaction m_transaction;  // the transaction providing isolation for the computation
    Iterator m_iterator; // an iterator attached to the transaction
    thread m_handle; // underlying thread

    // internal state for #process_vertex()
    uint64_t m_n1 =0; // current vertex being processed
    uint64_t m_n2 =0; // current neighbour being visited
    vector<uint64_t> m_neighbours; // neighbours of n1
    uint64_t m_marker =0; // current position in the neighbours vector, to merge shared neighbours
    uint64_t m_num_triangles =0; // current number of triangles found for `n1'

    // Process the given vertex
    void process_vertex(uint64_t vertex_id);

public:
    // Init
    LCC_Worker(LCC_Master* master, const Transaction& transaction);

    // Destructor
    ~LCC_Worker();

    // Main thread
    void execute();

    // Wait for the worker's thread to terminate
    void join();
};

/*****************************************************************************
 *                                                                           *
 *  LCC_Master                                                               *
 *                                                                           *
 *****************************************************************************/
#undef COUT_CLASS_NAME
#define COUT_CLASS_NAME "LCC_Master"

LCC_Master::LCC_Master(Teseo* teseo, Transaction transaction, utility::TimeoutService* timeout)
    : m_teseo(teseo), m_transaction(transaction), m_scores(nullptr), m_num_triangles(nullptr), m_next(0), m_timeout(timeout) {
}

LCC_Master::~LCC_Master() {
    delete[] m_scores; m_scores = nullptr;
    delete[] m_num_triangles; m_num_triangles = nullptr;
}

void LCC_Master::initialise(){
    assert(m_num_triangles == nullptr && "Already initialised");

    const uint64_t num_vertices = m_transaction.num_vertices();
    m_scores = new double[num_vertices](); // init to 0
    m_num_triangles = new atomic<uint64_t>[num_vertices](); // init to 0
}

void LCC_Master::compute_scores(){
    //common::Timer timer; timer.start();

    for(uint64_t i = 0, N = m_transaction.num_vertices(); i < N; i++){
        uint64_t num_triangles = m_num_triangles[i];
        if(num_triangles > 0){
            uint64_t degree = m_transaction.degree(i, /* logical ? */ true);
            uint64_t max_num_edges = degree * (degree -1);
            double score = static_cast<double>(num_triangles) / max_num_edges;
            m_scores[i] = score;
        } // else m_scores[i] = 0 (default value)
    }

    //timer.stop();
    //COUT_DEBUG_FORCE("compute_scores executed in: " << timer);
}

double* LCC_Master::execute() {
    common::Timer timer; timer.start();

    // init the state and the side information for each vertex
    initialise();
    //COUT_DEBUG_FORCE("Initialisation time: " << timer);

    // start the workers
    assert(LCC_NUM_WORKERS >= 1 && "At least one worker should be set");
    vector<LCC_Worker*> workers;
    workers.reserve(LCC_NUM_WORKERS);
    for(uint64_t worker_id = 0; worker_id < LCC_NUM_WORKERS; worker_id++ ){
        workers.push_back(new LCC_Worker(this, m_transaction));
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

bool LCC_Master::next_task(uint64_t* output_vtx_start /* inclusive */, uint64_t* output_vtx_end /* exclusive */) {
    uint64_t logical_start = m_next.fetch_add(LCC_TASK_SIZE); /* return the previous value of m_next */
    uint64_t num_vertices = m_transaction.num_vertices();
    if(logical_start >= num_vertices || m_timeout->is_timeout()){
        return false;
    } else {
        uint64_t logical_end = std::min(logical_start + LCC_TASK_SIZE, num_vertices);

        *output_vtx_start = logical_start;
        *output_vtx_end = logical_end;

        return true;
    }
}

Teseo* LCC_Master::teseo() {
    return m_teseo;
}

atomic<uint64_t>& LCC_Master::num_triangles(uint64_t vertex_id) {
    return m_num_triangles[vertex_id];
}

/*****************************************************************************
 *                                                                           *
 *  LCC_Worker                                                               *
 *                                                                           *
 *****************************************************************************/

#undef COUT_CLASS_NAME
#define COUT_CLASS_NAME "LCC_Worker"

LCC_Worker::LCC_Worker(LCC_Master* master, const Transaction& transaction) : m_master(master), m_transaction(transaction), m_iterator(m_transaction.iterator()) {
    m_handle = thread { &LCC_Worker::execute, this };
}

LCC_Worker::~LCC_Worker(){ }

void LCC_Worker::execute() {
    COUT_DEBUG("Worker started");

    m_master->teseo()->register_thread();

    uint64_t v_start, v_end;
    while(m_master->next_task(&v_start, &v_end)){
        for(uint64_t v = v_start; v < v_end; v++){
            process_vertex(v);
        }
    }

    m_iterator.close();
    m_master->teseo()->unregister_thread();

    COUT_DEBUG("Worker terminated");
}

void LCC_Worker::join(){
    m_handle.join();
}

void LCC_Worker::process_vertex(uint64_t n1) {
    COUT_DEBUG("vertex: " << n1);
    m_n1 = n1;
    m_neighbours.clear();
    m_num_triangles = 0; // reset the number of triangles

    m_iterator.edges(n1, /* logical ? */ true, [this](uint64_t n2){
        if(n2 > m_n1) return false; // we're done with n1
        m_n2 = n2;
        m_neighbours.push_back(n2);
        m_marker = 0; // reset the marker

        m_iterator.edges(n2, /* logical ? */ true, [this](uint64_t n3){
            if(n3 > m_n2) return false; // we're done with `n2'
            assert(m_n1 > m_n2 && m_n2 > n3); // we're looking for triangles of the kind c - b - a, with c > b && b > a

            COUT_DEBUG("  candidate: " << m_n1 << " - " << m_n2 << " - " << n3 << ", marker[" << m_marker << "] = " << m_neighbours[m_marker]);

            if(n3 > m_neighbours[m_marker]){ // merge with m_neighbours
                do {
                    m_marker ++;
                } while(m_marker < m_neighbours.size() && n3 > m_neighbours[m_marker]);
                if(m_marker >= m_neighbours.size()) return false; // there is nothing left to merge
            }

            if(n3 == m_neighbours[m_marker]){ // match !
                COUT_DEBUG("    match: " << m_n1 << " - " << m_n2 << " - " << n3 << " and " << m_n1 << " - " << n3 << " - " << m_n2);

                m_num_triangles += 2; // we've discovered both n1 - n2 - n3 and n1 - n3 - n2; with n1 > n2 > n3

                // increase the contribution for n2
                m_master->num_triangles(m_n2) += 2;

                // increase the contribution for n3
                m_master->num_triangles(n3) += 2;

                m_marker++;
                if(m_marker >= m_neighbours.size()) return false; // there is nothing left to merge
            }

            return true; // keep scanning

        });

        return true; // keep scanning
    });

    if(m_num_triangles != 0){
        m_master->num_triangles(m_n1) += m_num_triangles;
    }
}

} // anon namespace

#undef COUT_CLASS_NAME
#define COUT_CLASS_NAME "TeseoDriverLCC"

TeseoDriverLCC::TeseoDriverLCC(bool is_directed, bool read_only) : TeseoDriver(is_directed, read_only){

}

void TeseoDriverLCC::lcc(const char* dump2file){
    unique_ptr<double[]> scores;
    utility::TimeoutService timeout { m_timeout };
    common::Timer timer; timer.start();
    TESEO->register_thread();
    RegisterThread rt { TESEO } ;
    auto transaction = TESEO->start_transaction(/* read only ? */ m_read_only);

    { // restrict the scope to allow the dtor to clean up
        LCC_Master algorithm ( TESEO, transaction, &timeout );
        scores.reset( algorithm.execute() );
    }
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);  }

    // translate the vertex IDs
    uint64_t N = transaction.num_vertices();
    vector<pair<uint64_t, double>> external_ids(N);
    for (uint64_t v = 0; v < N; v++) {
        external_ids[v] = make_pair(transaction.vertex_id(v), scores[v]);
    }

    // store the results in the given file
    if(dump2file != nullptr)
        save_results(external_ids, dump2file);
}

} // namespace
