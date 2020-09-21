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

#include "stinger.hpp"

#include <cinttypes>
#include <fstream>
#include <limits>
#include <mutex>

#include "common/system.hpp"
#include "common/timer.hpp"
#include "third-party/gapbs/gapbs.hpp"
#include "utility/timeout_service.hpp"
#include "stinger_core/stinger.h"
#include "stinger_error.hpp"

using namespace gapbs;
using namespace libcuckoo;
using namespace std;

// Macros
#define STINGER reinterpret_cast<struct stinger*>(m_stinger_graph)
// Stinger stores the weights as int64_t, these macros handle the conversion to double
#define INT2DBL(v) ( *(reinterpret_cast<double*>(&(v))) )
#define DBL2INT(v) ( *(reinterpret_cast<int64_t*>(&(v))) )
#undef CURRENT_ERROR_TYPE
#define CURRENT_ERROR_TYPE ::gfe::library::StingerError


/******************************************************************************
 *                                                                            *
 *  Debug                                                                     *
 *                                                                            *
 *****************************************************************************/
//#define DEBUG
namespace gfe { extern mutex _log_mutex [[maybe_unused]]; }
#define COUT_DEBUG_FORCE(msg) { scoped_lock<mutex> lock(::gfe::_log_mutex); std::cout << "[Stinger::" << __FUNCTION__ << "] [" << common::concurrency::get_thread_id() << "] " << msg << std::endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif

namespace gfe::library {

/******************************************************************************
 *                                                                            *
 *  Init                                                                      *
 *                                                                            *
 *****************************************************************************/

StingerRef::StingerRef(bool directed): Stinger(directed) { }

StingerRef::~StingerRef() { }

/******************************************************************************
 *                                                                            *
 *  Helpers                                                                   *
 *                                                                            *
 *****************************************************************************/

template<typename T>
static void save0(cuckoohash_map<uint64_t, T>& result, const char* dump2file){
    assert(dump2file != nullptr);
    COUT_DEBUG("save the results to: " << dump2file)

    fstream handle(dump2file, ios_base::out);
    if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

    auto list_entries = result.lock_table();

    for(const auto& p : list_entries){
        handle << p.first << " " << p.second << "\n";
    }

    list_entries.unlock();
    handle.close();
}

/******************************************************************************
 *                                                                            *
 *  BFS                                                                       *
 *                                                                            *
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
int64_t BUStep(stinger_t* g, int64_t NV, pvector<int64_t>& distances, int64_t distance, Bitmap &front, Bitmap &next) {
    int64_t awake_count = 0;
    next.reset();
    #pragma omp parallel for reduction(+ : awake_count) schedule(dynamic, 1024)
    for (int64_t u=0; u < NV; u++) {
        if (distances[u] < 0){ // the node has not been visited yet
            bool skip = false;
            STINGER_FORALL_IN_EDGES_OF_VTX_BEGIN(g, u) { // check whether any of its incoming neighbours has already been visited
                if (!skip && front.get_bit(STINGER_EDGE_DEST)) {
//                    parent[u] = STINGER_EDGE_DEST; // if we were interested in the parent
                    distances[u] = distance; // on each BUStep, all nodes will have the same distance
                    awake_count++;
                    next.set_bit(u);
                    skip = true; // done, the node can be reached by the current frontier
                }
            } STINGER_FORALL_EDGES_OF_VTX_END();
        }
    }
    return awake_count;
}

static
int64_t TDStep(stinger_t* g, pvector<int64_t>& distances, int64_t distance, SlidingQueue<int64_t>& queue) {
    int64_t scout_count = 0;
    #pragma omp parallel
    {
        QueueBuffer<int64_t> lqueue(queue);
        #pragma omp for reduction(+ : scout_count) schedule(dynamic, 64)
        for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
            int64_t u = *q_iter;
            STINGER_FORALL_OUT_EDGES_OF_VTX_BEGIN(g, u) {
                int64_t curr_val = distances[STINGER_EDGE_DEST];
                if (curr_val < 0 && compare_and_swap(distances[STINGER_EDGE_DEST], curr_val, distance)) {
                    lqueue.push_back(STINGER_EDGE_DEST);
                    scout_count += -curr_val;
                }
            } STINGER_FORALL_OUT_EDGES_OF_VTX_END();
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
void BitmapToQueue(stinger_t* g, int64_t NV, const Bitmap &bm, SlidingQueue<int64_t> &queue) {
    #pragma omp parallel
    {
        QueueBuffer<int64_t> lqueue(queue);
        #pragma omp for
        for (int64_t n=0; n < NV; n++)
            if (bm.get_bit(n))
                lqueue.push_back(n);
        lqueue.flush();
    }
    queue.slide_window();
}

static
pvector<int64_t> InitDistances(stinger_t* g, int64_t NV){
    pvector<int64_t> distances(NV);
    #pragma omp parallel for
    for (uint64_t n = 0; n < NV; n++){
        int64_t out_degree = stinger_outdegree_get(g, n);
        distances[n] = out_degree != 0 ? - out_degree : -1;
    }
    return distances;
}

static
pvector<int64_t> do_bfs(stinger_t* g, int64_t NV, uint64_t num_out_edges, int64_t source, utility::TimeoutService& timer, int alpha = 15, int beta = 18) {
    // The implementation from GAP BS reports the parent (which indeed it should make more sense), while the one required by
    // Graphalytics only returns the distance
    pvector<int64_t> distances = InitDistances(g, NV);
    distances[source] = 0;

    SlidingQueue<int64_t> queue(NV);
    queue.push_back(source);
    queue.slide_window();
    Bitmap curr(NV);
    curr.reset();
    Bitmap front(NV);
    front.reset();
    int64_t edges_to_check = num_out_edges; //g.num_edges_directed();
    int64_t scout_count = stinger_outdegree_get(g, source);
    int64_t distance = 1; // current distance
    while (!timer.is_timeout() && !queue.empty()) {

        if (scout_count > edges_to_check / alpha) {
            int64_t awake_count, old_awake_count;
            QueueToBitmap(queue, front);
            awake_count = queue.size();
            queue.slide_window();
            do {
                old_awake_count = awake_count;
                awake_count = BUStep(g, NV, distances, distance, front, curr);
                front.swap(curr);
                distance++;
            } while ((awake_count >= old_awake_count) || (awake_count > NV / beta));
            BitmapToQueue(g, NV, front, queue);
            scout_count = 1;
        } else {
            edges_to_check -= scout_count;
            scout_count = TDStep(g, distances, distance, queue);
            queue.slide_window();
            distance++;
        }
    }

    return distances;
}

} // anon namespace

static void save_bfs(cuckoohash_map<uint64_t, int64_t>& result, const char* dump2file){
    assert(dump2file != nullptr);
    COUT_DEBUG("save the results to: " << dump2file)

    fstream handle(dump2file, ios_base::out);
    if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

    auto list_entries = result.lock_table();

    for(const auto& p : list_entries){
        handle << p.first << " ";

        // if  the vertex was not reached, the algorithm sets its distance to < 0
        if(p.second < 0){
            handle << numeric_limits<int64_t>::max();
        } else {
            handle << p.second;
        }
        handle << "\n";
    }

    list_entries.unlock();
    handle.close();
}

void StingerRef::bfs(uint64_t source_external_id, const char* dump2file){
    utility::TimeoutService tcheck ( chrono::seconds{ m_timeout } );
    common::Timer timer; timer.start();
//    uint64_t nv = stinger_max_active_vertex(STINGER); // it doesn't register the coefficient if the last vertices are isolated
    uint64_t nv = get_max_num_mappings(); // number of mappings
    int64_t source_internal_id = get_internal_id(source_external_id);

    // execute the algorithm
    auto result = do_bfs(STINGER, nv, num_edges(), source_internal_id, tcheck);
    if(tcheck.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // translate the vertex ids and set the distance of the non reached vertices to max()
    cuckoohash_map<uint64_t, int64_t> external_ids;
    #pragma omp parallel for
    for(uint64_t internal_id = 0; internal_id < nv; internal_id++){
        if(stinger_vtype_get(STINGER, internal_id) == 0){ // if = 1, the node is marked for deletion
            char* vertex_id_name = nullptr; uint64_t vertex_id_name_length = 0;
            int rc = stinger_mapping_physid_get(STINGER, internal_id, &vertex_id_name, &vertex_id_name_length);
            if( rc == 0 ){ // mapping found
                uint64_t external_id = stoull(vertex_id_name);
                COUT_DEBUG("external_id: " << external_id << ", internal_id: " << internal_id);

                external_ids.insert(external_id, result[internal_id]);
            }
            free(vertex_id_name); vertex_id_name = nullptr;
        }
    }
    if(tcheck.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    if(dump2file != nullptr)
        save_bfs(external_ids, dump2file);
}

/******************************************************************************
 *                                                                            *
 *  PageRank                                                                  *
 *                                                                            *
 *****************************************************************************/

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
pvector<double> do_pagerank(stinger_t* G, uint64_t num_active_vertices, uint64_t num_total_vertices, uint64_t num_iterations, double damping_factor, utility::TimeoutService& timer) {
    // init
    const double init_score = 1.0 / num_active_vertices;
    const double base_score = (1.0 - damping_factor) / num_active_vertices;
    pvector<double> scores(num_total_vertices, init_score);
    pvector<double> outgoing_contrib(num_total_vertices, 0.0);

    // pagerank iterations
    for(uint64_t iteration = 0; iteration < num_iterations && !timer.is_timeout(); iteration++){
        double dangling_sum = 0.0;

        // for each node, precompute its contribution to all of its outgoing neighbours and, if it's a sink,
        // add its rank to the `dangling sum' (to be added to all nodes).
        #pragma omp parallel for reduction(+:dangling_sum)
        for(uint64_t v = 0; v < num_total_vertices; v++){
            if(stinger_vtype_get(G, v) == 1) continue; // vertex marked for deletion
            uint64_t out_degree = stinger_outdegree(G, v);
            if(out_degree == 0){ // this is a sink
//                COUT_DEBUG("[" << iteration << "] dsum " << v << ": "<< dangling_sum << " + " << scores[v] << " -> " << (dangling_sum + scores[v]));
                dangling_sum += scores[v];
            } else {
                outgoing_contrib[v] = scores[v] / out_degree;
            }
        }

        dangling_sum /= num_active_vertices;
//        COUT_DEBUG("[" << iteration << "] base score: " << base_score << ", dangling sum: " << dangling_sum);

        // compute the new score for each node in the graph
        #pragma omp parallel for schedule(dynamic, 64)
        for(uint64_t v = 0; v < num_total_vertices; v++){
            if(stinger_vtype_get(G, v) == 1) continue; // vertex marked for deletion
            double incoming_total = 0;

            // compute the score for the current iteration looking at the incoming edges
            STINGER_FORALL_IN_EDGES_OF_VTX_BEGIN(G, v) {
                incoming_total += outgoing_contrib[STINGER_EDGE_DEST];
            } STINGER_FORALL_IN_EDGES_OF_VTX_END();

            // update the score
//            COUT_DEBUG("[" << iteration << "] contrib[" << v << "] : " << incoming_total);
            scores[v] = base_score + damping_factor * (incoming_total + dangling_sum);
//            COUT_DEBUG("[" << iteration << "] score[" << v << "] " << scores[v]);
        }

    }

    return scores;
}

void StingerRef::pagerank(uint64_t num_iterations, double damping_factor, const char* dump2file){
    utility::TimeoutService tcheck ( chrono::seconds{ m_timeout } );
    common::Timer timer; timer.start();
    const size_t num_registered_vertices = stinger_max_active_vertex(STINGER) +1; // logical vertex
    const size_t num_active_vertices = num_vertices();

    // execute the algorithm
    auto result = do_pagerank(STINGER, num_active_vertices, num_registered_vertices, num_iterations, damping_factor, tcheck);
    if(tcheck.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // translate the vertex ids and set the distance of the non reached vertices to max()
    cuckoohash_map<uint64_t, double> external_ids;
    #pragma omp parallel for
    for(uint64_t internal_id = 0; internal_id < num_registered_vertices; internal_id++){
        if(stinger_vtype_get(STINGER, internal_id) == 0){ // if = 1, the node is marked for deletion
            char* vertex_id_name = nullptr; uint64_t vertex_id_name_length = 0;
            int rc = stinger_mapping_physid_get(STINGER, internal_id, &vertex_id_name, &vertex_id_name_length);
            if( rc == 0 ){ // mapping found
                uint64_t external_id = stoull(vertex_id_name);
                COUT_DEBUG("external_id: " << external_id << ", internal_id: " << internal_id);

                external_ids.insert(external_id, result[internal_id]);
            }
            free(vertex_id_name); vertex_id_name = nullptr;
        }
    }
    if(tcheck.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    if(dump2file != nullptr)
        save0(external_ids, dump2file);
}

/******************************************************************************
 *                                                                            *
 *  Weakly connected components                                               *
 *                                                                            *
 *****************************************************************************/

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
pvector<uint64_t> do_wcc(stinger_t* G, uint64_t num_total_vertices, utility::TimeoutService& timer) {
    // init
    pvector<uint64_t> comp(num_total_vertices);
    #pragma omp parallel for
    for (uint64_t n = 0; n < num_total_vertices; n++){
        comp[n] = n;
    }

    bool change = true;
    while (change && !timer.is_timeout()) {
        change = false;

        #pragma omp parallel for schedule(dynamic, 64)
        for (uint64_t u = 0; u < num_total_vertices; u++){
            STINGER_FORALL_OUT_EDGES_OF_VTX_BEGIN(G, u) {
                uint64_t v = STINGER_EDGE_DEST;
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
            } STINGER_FORALL_OUT_EDGES_OF_VTX_END();
        }

        #pragma omp parallel for schedule(dynamic, 64)
        for (uint64_t n = 0; n < num_total_vertices; n++){
            while (comp[n] != comp[comp[n]]) {
                comp[n] = comp[comp[n]];
            }
        }
    }

    return comp;
}

void StingerRef::wcc(const char* dump2file) {
    // init
    utility::TimeoutService tcheck ( chrono::seconds{ m_timeout } );
    common::Timer timer; timer.start();
//    uint64_t nv = stinger_max_active_vertex(STINGER); // it doesn't register the coefficient if the last vertices are isolated
    uint64_t nv = get_max_num_mappings(); // number of mappings

    // execute the algorithm
    auto result = do_wcc(STINGER, nv, tcheck);
    if(tcheck.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // translate the vertex ids and set the distance of the non reached vertices to max()
    cuckoohash_map<uint64_t, uint64_t> external_ids;
    #pragma omp parallel for
    for(uint64_t internal_id = 0; internal_id < nv; internal_id++){
        if(stinger_vtype_get(STINGER, internal_id) == 0){ // if = 1, the node is marked for deletion
            char* vertex_id_name = nullptr; uint64_t vertex_id_name_length = 0;
            int rc = stinger_mapping_physid_get(STINGER, internal_id, &vertex_id_name, &vertex_id_name_length);
            if( rc == 0 ){ // mapping found
                uint64_t external_id = stoull(vertex_id_name);
                COUT_DEBUG("external_id: " << external_id << ", internal_id: " << internal_id);

                external_ids.insert(external_id, result[internal_id]);
            }
            free(vertex_id_name); vertex_id_name = nullptr;
        }
    }
    if(tcheck.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    if(dump2file != nullptr)
        save0(external_ids, dump2file);
}


/******************************************************************************
 *                                                                            *
 *  Single Source Shortest Path (SSSP)                                        *
 *                                                                            *
 *****************************************************************************/

static const size_t kMaxBin = numeric_limits<size_t>::max()/2;

static pvector<double> do_sssp(stinger_t* G, uint64_t num_total_vertices, uint64_t num_total_edges, uint64_t source, double delta, utility::TimeoutService& timer){
    // Init
    pvector<double> dist(num_total_vertices, numeric_limits<double>::infinity());
    dist[source] = 0;
    pvector<uint64_t> frontier(num_total_edges);
    // two element arrays for double buffering curr=iter&1, next=(iter+1)&1
    size_t shared_indexes[2] = {0, kMaxBin};
    size_t frontier_tails[2] = {1, 0};
    frontier[0] = source;

    #pragma omp parallel
    {
        vector<vector<uint64_t>> local_bins(0);
        size_t iter = 0;

        while (shared_indexes[iter&1] != kMaxBin) {
            size_t &curr_bin_index = shared_indexes[iter&1];
            size_t &next_bin_index = shared_indexes[(iter+1)&1];
            size_t &curr_frontier_tail = frontier_tails[iter&1];
            size_t &next_frontier_tail = frontier_tails[(iter+1)&1];
            #pragma omp for nowait schedule(dynamic, 64)
            for (size_t i = 0; i < curr_frontier_tail; i++) {
                uint64_t u = frontier[i];
                COUT_DEBUG("[" << iter << "] examine " << u);
                if (dist[u] >= delta * static_cast<double>(curr_bin_index)) {

                    STINGER_FORALL_OUT_EDGES_OF_VTX_BEGIN(G, u) {
                        uint64_t v = STINGER_EDGE_DEST;
                        double w = INT2DBL(STINGER_EDGE_WEIGHT);

                        double old_dist = dist[v];
                        double new_dist = dist[u] + w;
                        if(new_dist < old_dist){
                            bool changed_dist = true;
                            while (!compare_and_swap(dist[v], old_dist, new_dist)) {
                                old_dist = dist[v];
                                if (old_dist <= new_dist) {
                                    changed_dist = false;
                                    break;
                                }
                            }

                            if(changed_dist){
                                size_t dest_bin = new_dist/delta;
                                COUT_DEBUG("Update " << v << " from " << old_dist << " to " << new_dist);
                                if (dest_bin >= local_bins.size()) {
                                    local_bins.resize(dest_bin+1);
                                }
                                local_bins[dest_bin].push_back(v);
                            }
                        }
                    } STINGER_FORALL_OUT_EDGES_OF_VTX_END();

                }
            }

            for (size_t i = curr_bin_index; i < local_bins.size(); i++) {
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
                size_t copy_start = fetch_and_add(next_frontier_tail, local_bins[next_bin_index].size());
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

void StingerRef::sssp(uint64_t source_external_id, const char* dump2file){
    utility::TimeoutService tcheck ( chrono::seconds{ m_timeout } );
    common::Timer timer; timer.start();
//    uint64_t nv = stinger_max_active_vertex(STINGER); // it doesn't register the coefficient if the last vertices are isolated
    uint64_t nv = get_max_num_mappings(); // number of mappings
    uint64_t ne = num_edges();
    int64_t source_internal_id = get_internal_id(source_external_id);

    // execute the algorithm
    double delta = 2.0; // same value used in the GAPBS, at least for most graphs
    auto result = do_sssp(STINGER, nv, ne, source_internal_id, delta, tcheck);
    if(tcheck.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // translate the vertex ids and set the distance of the non reached vertices to max()
    cuckoohash_map<uint64_t, double> external_ids;
    #pragma omp parallel for
    for(uint64_t internal_id = 0; internal_id < nv; internal_id++){
        if(stinger_vtype_get(STINGER, internal_id) == 0){ // if = 1, the node is marked for deletion
            char* vertex_id_name = nullptr; uint64_t vertex_id_name_length = 0;
            int rc = stinger_mapping_physid_get(STINGER, internal_id, &vertex_id_name, &vertex_id_name_length);
            if( rc == 0 ){ // mapping found
                uint64_t external_id = stoull(vertex_id_name);
                COUT_DEBUG("external_id: " << external_id << ", internal_id: " << internal_id);

                external_ids.insert(external_id, result[internal_id]);
            }
            free(vertex_id_name); vertex_id_name = nullptr;
        }
    }
    if(tcheck.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    if(dump2file != nullptr)
        save0(external_ids, dump2file);
}


} // namespace
