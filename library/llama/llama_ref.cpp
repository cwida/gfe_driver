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

#include "llama_ref.hpp"

#include <cinttypes>
#include <fstream>
#include <limits>
#include <mutex>

#include "common/system.hpp"
#include "common/timer.hpp"
#include "third-party/gapbs/gapbs.hpp"
#include "utility/timeout_service.hpp"
#include "llama_internal.hpp"

using namespace gapbs;
using namespace std;

/******************************************************************************
 *                                                                            *
 *  Debug                                                                     *
 *                                                                            *
 *****************************************************************************/
//#define DEBUG
namespace gfe { extern mutex _log_mutex [[maybe_unused]]; }
#define COUT_DEBUG_FORCE(msg) { scoped_lock<mutex> lock(::gfe::_log_mutex); std::cout << "[LLAMARef::" << __FUNCTION__ << "] [" << common::concurrency::get_thread_id() << "] " << msg << std::endl; }
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

LLAMARef::LLAMARef(bool directed): LLAMAClass(directed) { }

LLAMARef::~LLAMARef() { }

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


static
uint64_t _llama_outdegree(ll_mlcsr_ro_graph& snapshot, bool is_directed, int64_t llama_vertex_id){
    if(is_directed){
        return snapshot.out_degree(llama_vertex_id);
    } else {
        return snapshot.out_degree(llama_vertex_id) + snapshot.in_degree(llama_vertex_id);
    }
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
int64_t BUStep(ll_mlcsr_ro_graph& graph, bool is_directed, pvector<int64_t>& distances, int64_t distance, Bitmap &front, Bitmap &next) {
    int64_t awake_count = 0;
    next.reset();
    #pragma omp parallel for reduction(+ : awake_count) schedule(dynamic, 1024)
    for (int64_t u=0; u < graph.max_nodes(); u++) {
        if (distances[u] < 0){ // the node has not been visited yet

            ll_edge_iterator iterator;
            graph.in_iter_begin_fast(iterator, u);
            bool done = false;
            for (edge_t e = graph.in_iter_next_fast(iterator); e != LL_NIL_EDGE && !done; e = graph.in_iter_next_fast(iterator)) {
                node_t n = LL_ITER_OUT_NEXT_NODE(*G, iterator, e);

                if (front.get_bit(n)) {
                    distances[u] = distance; // on each BUStep, all nodes will have the same distance
                    awake_count++;
                    next.set_bit(u);
                    done = true;
                }
            }

            if(!is_directed && !done){ // the actual set of edges in undirected graphs is composed by both LLAMA's incoming and outgoing edges
                graph.out_iter_begin(iterator, u);
                for (edge_t e = graph.out_iter_next(iterator); e != LL_NIL_EDGE && !done; e = graph.out_iter_next(iterator)) {
                    node_t n = LL_ITER_OUT_NEXT_NODE(*G, iterator, e);
                    if (front.get_bit(n)) {
                        distances[u] = distance; // on each BUStep, all nodes will have the same distance
                        awake_count++;
                        next.set_bit(u);
                        done = true;
                    }
                }
            }
        }
    }
    return awake_count;
}

static
int64_t TDStep(ll_mlcsr_ro_graph& graph, bool is_directed, pvector<int64_t>& distances, int64_t distance, SlidingQueue<int64_t>& queue) {
    int64_t scout_count = 0;
    #pragma omp parallel
    {
        QueueBuffer<int64_t> lqueue(queue);
        #pragma omp for reduction(+ : scout_count)
        for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
            int64_t u = *q_iter;

            ll_edge_iterator iterator;
            graph.out_iter_begin(iterator, u);
            for (edge_t e = graph.out_iter_next(iterator); e != LL_NIL_EDGE; e = graph.out_iter_next(iterator)) {
                node_t n = LL_ITER_OUT_NEXT_NODE(*G, iterator, e);

                int64_t curr_val = distances[n];
                if (curr_val < 0 && compare_and_swap(distances[n], curr_val, distance)) {
                    lqueue.push_back(n);
                    scout_count += -curr_val;
                }
            }

            if(!is_directed){ // the actual set of edges in undirected graphs is composed by both LLAMA's incoming and outgoing edges
                graph.in_iter_begin_fast(iterator, u);
                for (edge_t e = graph.in_iter_next_fast(iterator); e != LL_NIL_EDGE; e = graph.in_iter_next_fast(iterator)) {
                    node_t n = LL_ITER_OUT_NEXT_NODE(*G, iterator, e);

                    int64_t curr_val = distances[n];
                    if (curr_val < 0 && compare_and_swap(distances[n], curr_val, distance)) {
                        lqueue.push_back(n);
                        scout_count += -curr_val;
                    }
                }
            }
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
void BitmapToQueue(ll_mlcsr_ro_graph& graph, const Bitmap &bm, SlidingQueue<int64_t> &queue) {
    #pragma omp parallel
    {
        QueueBuffer<int64_t> lqueue(queue);
        #pragma omp for
        for (int64_t n=0; n < graph.max_nodes(); n++)
            if (bm.get_bit(n))
                lqueue.push_back(n);
        lqueue.flush();
    }
    queue.slide_window();
}

static
pvector<int64_t> InitDistances(ll_mlcsr_ro_graph& graph, bool is_directed){
    pvector<int64_t> distances(graph.max_nodes());
    #pragma omp parallel for
    for (int64_t n = 0; n < graph.max_nodes(); n++){
        int64_t out_degree = _llama_outdegree(graph, is_directed, n);
        distances[n] = out_degree != 0 ? - out_degree : -1;
    }
    return distances;
}

static
pvector<int64_t> do_bfs(ll_mlcsr_ro_graph& graph, bool is_directed, uint64_t num_edges, int64_t source, utility::TimeoutService& timer, int alpha = 15, int beta = 18) {
    // The implementation from GAP BS reports the parent (which indeed it should make more sense), while the one required by
    // Graphalytics only returns the distance
    pvector<int64_t> distances = InitDistances(graph, is_directed);
    distances[source] = 0;

    SlidingQueue<int64_t> queue(graph.max_nodes());
    queue.push_back(source);
    queue.slide_window();
    Bitmap curr(graph.max_nodes());
    curr.reset();
    Bitmap front(graph.max_nodes());
    front.reset();
    int64_t edges_to_check = num_edges; //g.num_edges_directed();
    int64_t scout_count = _llama_outdegree(graph, is_directed, source);
    int64_t distance = 1; // current distance
    while (!timer.is_timeout() && !queue.empty()) {

        if (scout_count > edges_to_check / alpha) {
            int64_t awake_count, old_awake_count;
            QueueToBitmap(queue, front);
            awake_count = queue.size();
            queue.slide_window();
            do {
                old_awake_count = awake_count;
                awake_count = BUStep(graph, is_directed, distances, distance, front, curr);
                front.swap(curr);
                distance++;
            } while ((awake_count >= old_awake_count) || (awake_count > graph.max_nodes() / beta));
            BitmapToQueue(graph, front, queue);
            scout_count = 1;
        } else {
            edges_to_check -= scout_count;
            scout_count = TDStep(graph, is_directed, distances, distance, queue);
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

void LLAMARef::bfs(uint64_t external_source_vertex_id, const char* dump2file){
    utility::TimeoutService tcheck { m_timeout };
    common::Timer timer; timer.start();

    // retrieve the latest snapshot the internal source_vertex_id
    shared_lock<shared_mutex> slock(m_lock_checkpoint);
    auto graph = get_snapshot();
    int64_t llama_source_vertex_id;
    if (! m_vmap_read_only.find(external_source_vertex_id, llama_source_vertex_id) ){ // side effect, it assigns llama_source_vertex_id
        ERROR("The given vertex does not exist: " << external_source_vertex_id);
    }
    uint64_t num_edges = 0;
    for(uint64_t i = 0; i < graph.num_levels(); i++){ num_edges += graph.max_edges(i); }

    slock.unlock(); // here we lose the ability to refer to m_vmap_read_only from now on

    // execute the BFS algorithm
    auto result = do_bfs(graph, is_directed(), num_edges, llama_source_vertex_id, tcheck);
    if(tcheck.is_timeout()){  RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // translate from llama vertex ids to external vertex ids
    auto names = graph.get_node_property_64(g_llama_property_names);
    assert(names != nullptr && "Wrong string ID to refer the property attached to the vertices");
    cuckoohash_map</* external id */ uint64_t, /* distance */ int64_t> external_ids;
    #pragma omp parallel for
    for(node_t llama_node_id = 0; llama_node_id < graph.max_nodes(); llama_node_id++){
        // first, does this node exist (or it's a gap?)
        // this is a bit of a stretch: the impl~ from llama assumes that a node does not exist only if it does not have any incoming or outgoing edges.
        if(!graph.node_exists(llama_node_id)) continue;

        // second, what's it's real node ID, in the external domain (e.g. user id)
        uint64_t external_node_id = names->get(llama_node_id);

        // third, its distance
        auto distance = result[llama_node_id];

        // finally, register the association
        external_ids.insert(external_node_id, distance);
    }

    if(tcheck.is_timeout()) RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);

    // store the results in the given file
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
pvector<double> do_pagerank(ll_mlcsr_ro_graph& graph, uint64_t current_num_vertices, bool is_directed, uint64_t num_iterations, double damping_factor, utility::TimeoutService& timer) {
    // init
    const double init_score = 1.0 / current_num_vertices;
    const double base_score = (1.0 - damping_factor) / current_num_vertices;
    pvector<double> scores(graph.max_nodes(), init_score);
    pvector<double> outgoing_contrib(graph.max_nodes(), 0.0);

    // pagerank iterations
    for(uint64_t iteration = 0; iteration < num_iterations && !timer.is_timeout(); iteration++){
        double dangling_sum = 0.0;

        // for each node, precompute its contribution to all of its outgoing neighbours and, if it's a sink,
        // add its rank to the `dangling sum' (to be added to all nodes).
        #pragma omp parallel for reduction(+:dangling_sum)
        for(uint64_t v = 0; v < graph.max_nodes(); v++){
            // okay, LLAMA does not know if a vertex does not exist or it's just a gap / empty cell in its arrays
            // we assume all nodes with no incoming or outgoing edges are gaps
            if(graph.out_degree(v) + graph.in_degree(v) == 0) continue;

            uint64_t out_degree = _llama_outdegree(graph, is_directed, v);
            if(out_degree == 0){ // this is a sink
                dangling_sum += scores[v];
            } else {
                outgoing_contrib[v] = scores[v] / out_degree;
            }
        }

        dangling_sum /= current_num_vertices;

        // compute the new score for each node in the graph
        #pragma omp parallel for schedule(dynamic, 64)
        for(uint64_t v = 0; v < graph.max_nodes(); v++){
            // okay, LLAMA does not know if a vertex does not exist or it's just a gap / empty cell in its arrays
            // we assume all nodes with no incoming or outgoing edges are gaps
            if(graph.out_degree(v) + graph.in_degree(v) == 0) continue;
            double incoming_total = 0;

            // compute the score for the current iteration looking at the incoming edges
            ll_edge_iterator iterator;
            graph.in_iter_begin_fast(iterator, v);
            for (edge_t e = graph.in_iter_next_fast(iterator); e != LL_NIL_EDGE; e = graph.in_iter_next_fast(iterator)) {
                node_t n = LL_ITER_OUT_NEXT_NODE(*G, iterator, e);
                incoming_total += outgoing_contrib[n];
            }

            if(!is_directed){ // the actual set of edges in undirected graphs is composed by both LLAMA's incoming and outgoing edges
                graph.out_iter_begin(iterator, v);
                for (edge_t e = graph.out_iter_next(iterator); e != LL_NIL_EDGE; e = graph.out_iter_next(iterator)) {
                    node_t n = LL_ITER_OUT_NEXT_NODE(*G, iterator, e);
                    incoming_total += outgoing_contrib[n];
                }
            }

            // update the score
            scores[v] = base_score + damping_factor * (incoming_total + dangling_sum);
        }

    }

    return scores;
}

void LLAMARef::pagerank(uint64_t num_iterations, double damping_factor, const char* dump2file){
    utility::TimeoutService tcheck { m_timeout };
    common::Timer timer; timer.start();

    // retrieve the latest snapshot
    shared_lock<shared_mutex> slock(m_lock_checkpoint);
    auto graph = get_snapshot();
//    dump_snapshot(graph);
    auto current_num_vertices = num_vertices();
    slock.unlock();

    // execute the algorithm
    auto result = do_pagerank(graph, current_num_vertices, is_directed(), num_iterations, damping_factor, tcheck);
    if(tcheck.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // retrieve the external node ids
    auto names = graph.get_node_property_64(g_llama_property_names);
    assert(names != nullptr && "Wrong string ID to refer the property attached to the vertices");
    cuckoohash_map</* external id */ uint64_t, /* score */ double> external_ids;
    #pragma omp parallel for
    for(node_t llama_node_id = 0; llama_node_id < graph.max_nodes(); llama_node_id++){
        // first, does this node exist (or it's a gap?)
        // this is a bit of a stretch: the impl~ from llama assumes that a node does not exist only if it does not have any incoming or outgoing edges.
        if(!graph.node_exists(llama_node_id)) continue;

        // second, what's it's real node ID, in the external domain (e.g. user id)
        uint64_t external_id = names->get(llama_node_id);

        external_ids.insert(external_id, result[llama_node_id]);
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
pvector<uint64_t> do_wcc(ll_mlcsr_ro_graph& graph, bool is_directed, utility::TimeoutService& timer) {
    // init
    pvector<uint64_t> comp(graph.max_nodes());
    #pragma omp parallel for
    for (node_t n = 0; n < graph.max_nodes(); n++){
        comp[n] = n;
    }

    bool change = true;
    while (change && !timer.is_timeout()) {
        change = false;

        #pragma omp parallel for
        for (node_t u = 0; u < graph.max_nodes(); u++){
            ll_edge_iterator iterator;
            graph.out_iter_begin(iterator, u);
            for (edge_t e = graph.out_iter_next(iterator); e != LL_NIL_EDGE; e = graph.out_iter_next(iterator)) {
                node_t v = LL_ITER_OUT_NEXT_NODE(*G, iterator, e);
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

            if(!is_directed){ // the actual set of edges in undirected graphs is composed by both LLAMA's incoming and outgoing edges
                graph.in_iter_begin_fast(iterator, u);
                for (edge_t e = graph.in_iter_next_fast(iterator); e != LL_NIL_EDGE; e = graph.in_iter_next_fast(iterator)) {
                    node_t v = LL_ITER_OUT_NEXT_NODE(*G, iterator, e);
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
        }

        #pragma omp parallel for
        for (node_t n = 0; n < graph.max_nodes(); n++){
            while (comp[n] != comp[comp[n]]) {
                comp[n] = comp[comp[n]];
            }
        }
    }

    return comp;
}

void LLAMARef::wcc(const char* dump2file) {
    // init
    utility::TimeoutService tcheck { m_timeout };
    common::Timer timer; timer.start();

    // retrieve the latest snapshot the internal source_vertex_id
    shared_lock<shared_mutex> slock(m_lock_checkpoint);
    auto graph = get_snapshot();
    slock.unlock();

    // execute the algorithm
    auto result = do_wcc(graph, is_directed(), tcheck);
    if(tcheck.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // retrieve the external node ids
    auto names = graph.get_node_property_64(g_llama_property_names);
    assert(names != nullptr && "Wrong string ID to refer the property attached to the vertices");
    cuckoohash_map</* external id */ uint64_t, /* component */ uint64_t> external_ids;
    #pragma omp parallel for
    for(node_t llama_node_id = 0; llama_node_id < graph.max_nodes(); llama_node_id++){
        // first, does this node exist (or it's a gap?)
        // this is a bit of a stretch: the impl~ from llama assumes that a node does not exist only if it does not have any incoming or outgoing edges.
        if(!graph.node_exists(llama_node_id)) continue;

        // second, what's it's real node ID, in the external domain (e.g. user id)
        uint64_t external_id = names->get(llama_node_id);

        external_ids.insert(external_id, result[llama_node_id]);
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

static pvector<double> do_sssp(ll_mlcsr_ro_graph& graph, bool is_directed, uint64_t source, double delta, utility::TimeoutService& timer){
    // Init
    uint64_t num_total_vertices = graph.max_nodes();
    uint64_t num_total_edges = 0;
    for(int i = 0; i < graph.num_levels(); i++){ num_total_edges += graph.max_edges(i); }
    ll_mlcsr_edge_property<double>& weights = *reinterpret_cast<ll_mlcsr_edge_property<double>*>(graph.get_edge_property_64(g_llama_property_weights));
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

                    ll_edge_iterator iterator;
                    graph.out_iter_begin(iterator, u);
                    for (edge_t e = graph.out_iter_next(iterator); e != LL_NIL_EDGE; e = graph.out_iter_next(iterator)) {
                        node_t v = LL_ITER_OUT_NEXT_NODE(*G, iterator, e);
                        double w = weights[e];
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
                                COUT_DEBUG("Update " << v << " from " << old_dist << " to " << new_dist);
                                size_t dest_bin = new_dist/delta;
                                if (dest_bin >= local_bins.size()) {
                                    local_bins.resize(dest_bin+1);
                                }
                                local_bins[dest_bin].push_back(v);
                            }
                        }
                    }

                    if(!is_directed){ // the actual set of edges in undirected graphs is composed by both LLAMA's incoming and outgoing edges
                        graph.in_iter_begin_fast(iterator, u);
                        for (edge_t s_idx = graph.in_iter_next_fast(iterator); s_idx != LL_NIL_EDGE; s_idx = graph.in_iter_next_fast(iterator)) {
                            node_t v = LL_ITER_OUT_NEXT_NODE(*G, iterator, e);
                            edge_t e = graph.in_to_out(s_idx); // the weight is only a property of the outgoing edges
                            double w = weights[e];
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
                                    COUT_DEBUG("Update " << v << " from " << old_dist << " to " << new_dist);
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

void LLAMARef::sssp(uint64_t external_source_vertex_id, const char* dump2file){
    utility::TimeoutService tcheck { m_timeout };
    common::Timer timer; timer.start();

    // retrieve the latest snapshot the internal source_vertex_id
    shared_lock<shared_mutex> slock(m_lock_checkpoint);
    auto graph = get_snapshot();
//    dump_snapshot(graph);
    int64_t llama_source_vertex_id;
    if (! m_vmap_read_only.find(external_source_vertex_id, llama_source_vertex_id) ){ // side effect, it assigns llama_source_vertex_id
        ERROR("The given vertex does not exist: " << external_source_vertex_id);
    }
    slock.unlock();

    // execute the algorithm
    double delta = 2.0; // same value used in the GAPBS, at least for most graphs
    auto result = do_sssp(graph, is_directed(), llama_source_vertex_id, delta, tcheck);
    if(tcheck.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // translate the vertex ids
    auto names = graph.get_node_property_64(g_llama_property_names);
    assert(names != nullptr && "Wrong string ID to refer the property attached to the vertices");
    cuckoohash_map</* external id */ uint64_t, /* distance */ double> external_ids;
    #pragma omp parallel for
    for(node_t llama_node_id = 0; llama_node_id < graph.max_nodes(); llama_node_id++){
        // first, does this node exist (or it's a gap?)
        // this is a bit of a stretch: the impl~ from llama assumes that a node does not exist only if it does not have any incoming or outgoing edges.
        if(!graph.node_exists(llama_node_id)) continue;

        // second, what's it's real node ID, in the external domain (e.g. user id)
        uint64_t external_id = names->get(llama_node_id);

        COUT_DEBUG("Node map " << llama_node_id << " -> " << external_id);

        external_ids.insert(external_id, result[llama_node_id]);
    }
    if(tcheck.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    if(dump2file != nullptr)
        save0(external_ids, dump2file);
}


} // namespace


