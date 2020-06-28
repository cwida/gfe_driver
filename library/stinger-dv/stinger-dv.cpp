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

#include "stinger-dv.hpp"

#include <cassert>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <mutex>
#include <queue>

#include "stinger_core/stinger.h"

extern "C" {
#include "stinger_alg/weakly_connected_components.h"
}

#include "common/timer.hpp"
#include "library/stinger/stinger_error.hpp"
#include "utility/timeout_service.hpp"

using namespace std;

#define STINGER reinterpret_cast<struct stinger*>(m_stinger_graph)
#define INT2DBL(v) ( *(reinterpret_cast<double*>(&(v))) )
#define DBL2INT(v) ( *(reinterpret_cast<int64_t*>(&(v))) )

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG
extern mutex _log_mutex [[maybe_unused]];
#define COUT_DEBUG_FORCE(msg) { std::scoped_lock<std::mutex> lock{::_log_mutex}; std::cout << "[StingerDV::" << __FUNCTION__ << "] " << msg << std::endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif

/*****************************************************************************
 *                                                                           *
 *  Error                                                                    *
 *                                                                           *
 *****************************************************************************/
#undef CURRENT_ERROR_TYPE
#define CURRENT_ERROR_TYPE ::gfe::library::StingerError

namespace gfe::library {

/*****************************************************************************
 *                                                                           *
 *  Initialisation                                                           *
 *                                                                           *
 *****************************************************************************/

StingerDV::StingerDV(bool directed) : m_directed(directed){
    struct stinger_config_t config;
    memset(&config, 0, sizeof(config)); // init
    config.nv = 1ll<<32; // max number of vertices, 4G
    config.nebs = 0; // max number of edge blocks, 0=auto
    config.netypes = 1; // number of edge types, we are not going to use this feature, 1
    config.nvtypes = 1; // number of vertex types, we are not going to use this feature, 1
    config.no_map_none_etype = config.no_map_none_vtype = 1; // I think this is whether we want to attach names (as strings) to the vertex/edge types, 1 = we don't use this feature
    config.no_resize = 0; // when building stinger for the first time, allow to resize the internal structs if they do not fit into memory, 0 = allow resize

    m_stinger_graph = stinger_new_full(&config);
    assert(m_stinger_graph != nullptr && "Stinger allocation");
    if(m_stinger_graph == nullptr) ERROR("Cannot allocate stinger graph");
}

StingerDV::~StingerDV(){
    stinger_free_all(STINGER); m_stinger_graph = nullptr;
}

/******************************************************************************
 *                                                                            *
 *  Properties                                                                *
 *                                                                            *
 *****************************************************************************/
uint64_t StingerDV::num_edges() const {
    uint64_t num_directed_edges = static_cast<uint64_t>(stinger_total_edges(STINGER));
    if(!m_directed) num_directed_edges /= 2; /* divided by two, because we always store the edge for both source/destination */
    return num_directed_edges;
}

uint64_t StingerDV::num_vertices() const {
    return stinger_num_active_vertices(STINGER); // stinger_num_active_vertices only checks for vertices that have not been created
}

bool StingerDV::has_vertex(uint64_t vertex_id) const {
    return vertex_id < STINGER->max_nv;
}

double StingerDV::get_weight(uint64_t source, uint64_t destination) const {
    constexpr double NaN { numeric_limits<double>::signaling_NaN() };

    STINGER_FORALL_OUT_EDGES_OF_VTX_BEGIN(STINGER, source) {
        if(STINGER_EDGE_DEST == destination){
            return INT2DBL(STINGER_EDGE_WEIGHT);
        }
    } STINGER_FORALL_OUT_EDGES_OF_VTX_END();

    return NaN;
}

bool StingerDV::is_directed() const {
    return m_directed;
}

void StingerDV::set_timeout(uint64_t seconds) {
    m_timeout = seconds;
}

/******************************************************************************
 *                                                                            *
 *  Updates                                                                   *
 *                                                                            *
 *****************************************************************************/
bool StingerDV::add_vertex(uint64_t vertex_id){
    if(vertex_id >= STINGER->max_nv) ERROR("Invalid vertex ID: " << vertex_id << ", max vertex id allowed: " << STINGER->max_nv);
    return true; // lie
}

bool StingerDV::remove_vertex(uint64_t vertex_id){
    stinger_remove_vertex(STINGER, vertex_id); // ignore rc
    return true; // lie
}

bool StingerDV::add_edge(graph::WeightedEdge e){
    COUT_DEBUG("edge: " << e);

    int rc (0);
    if(m_directed){ // directed graph
        rc = stinger_insert_edge (STINGER, /* type, ignore */ 0 , e.source(), e.destination(), e.weight(), /* timestamp, ignore */ 0);
    } else { // undirected graph
        rc = stinger_insert_edge_pair(STINGER, /* type, ignore */ 0, e.source(), e.destination(), e.weight(), /* timestamp, ignore */ 0);
    }

    return rc >= 0;
}


bool StingerDV::add_edge_v2(gfe::graph::WeightedEdge e){
    add_vertex(e.source());
    add_vertex(e.destination());
    return add_edge(e);
}

bool StingerDV::remove_edge(graph::Edge e){
    int rc = 0;
    if(m_directed){ // directed graph
        // @return 1 on success, 0 if the edge is not found.
        rc = stinger_remove_edge(STINGER, /* type, ignore */ 0, e.source(), e.destination());

    } else { // undirected graph
        // @return < 0 in case of error, 0 edges not found, > 0 edges removed
        rc = stinger_remove_edge_pair(STINGER, /* type, ignore */ 0, e.source(), e.destination());
    }

    return rc > 0;
}

/******************************************************************************
 *                                                                            *
 *  Community detection through label propagation                             *
 *                                                                            *
 ******************************************************************************/
void StingerDV::cdlp(uint64_t max_iterations, const char* dump2file){
    auto tcheck = make_unique<utility::TimeoutService>( chrono::seconds{ m_timeout } );
    common::Timer timer; timer.start();

    bool change = true;
    uint64_t N = stinger_max_active_vertex(STINGER) +1;
    auto ptr_labels0 = make_unique<int64_t[]>(N); int64_t* labels0 = ptr_labels0.get();
    auto ptr_labels1 = make_unique<int64_t[]>(N); int64_t* labels1 = ptr_labels1.get();

    // initialisation
    #pragma omp parallel for
    for(int64_t n = 0; n < N; n++){
        labels0[n] = n;
    }

    // algorithm pass
    for(uint64_t i = 0; i < max_iterations && change; i++){
        COUT_DEBUG("iteration: " << i);
        if(tcheck->is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }
        change = false;

        #pragma omp parallel for shared(change)
        for(int64_t n = 0; n < N; n++){
            labels1[n] = cdlp_propagate(n, labels0);
            change |= (labels0[n] != labels1[n]);
        }

        std::swap(labels0, labels1); // next iteration
    }

    // save the computation into `dump2file'
    save(labels0, N, dump2file);
}

int64_t StingerDV::cdlp_propagate(int64_t vertex_id, int64_t* __restrict labels){
    unordered_map<int64_t, uint64_t> histogram;

    // compute the histogram
    STINGER_FORALL_OUT_EDGES_OF_VTX_BEGIN(STINGER, vertex_id) {
        histogram[labels[STINGER_EDGE_DEST]]++;
    } STINGER_FORALL_OUT_EDGES_OF_VTX_END();

    // cfr. Spec v0.9 pp 14 "If the graph is directed and a neighbor is reachable via both an incoming and
    // outgoing edge, its label will be counted twice"
    if(m_directed){
        STINGER_FORALL_IN_EDGES_OF_VTX_BEGIN(STINGER, vertex_id) {
            histogram[labels[STINGER_EDGE_DEST]]++;
        } STINGER_FORALL_IN_EDGES_OF_VTX_END();
    }

    // get the max label
    int64_t label_max = numeric_limits<int64_t>::max();
    uint64_t count_max = 0;
    for(const auto pair : histogram){
        if(pair.second > count_max || (pair.second == count_max && pair.first < label_max)){
            label_max = pair.first;
            count_max = pair.second;
        }
    }

    return label_max;
}

/******************************************************************************
 *                                                                            *
 *  Weakly connected components                                               *
 *                                                                            *
 *****************************************************************************/
void StingerDV::wcc(const char* dump2file) {
    // ignore the timeout as we use the impl~ from stinger
    int64_t nv = stinger_max_nv(STINGER);
    COUT_DEBUG("nv: " << nv);
    auto ptr_component_map = make_unique<int64_t[]>(nv); int64_t* component_map = ptr_component_map.get();
    parallel_shiloach_vishkin_components_of_type(STINGER, component_map, /* type, ignore */ 0); // already implemented in Stinger

    // store the final results (if required)
    save(component_map, stinger_max_active_vertex(STINGER) +1, dump2file);
}

/******************************************************************************
 *                                                                            *
 *  LCC                                                                       *
 *                                                                            *
 ******************************************************************************/
// Implementation based on stinger_alg/src/clustering.c

static int compare(const void* a, const void* b){
    return ( *(int64_t*)a - *(int64_t*)b );
}


void StingerDV::lcc(const char* dump2file){
    auto tcheck = make_unique<utility::TimeoutService>( chrono::seconds{ m_timeout } );
    common::Timer timer; timer.start();
//    uint64_t N = stinger_max_nv(STINGER); // assumes static vertices, this is 4G as currently set in the config
    uint64_t N = stinger_max_active_vertex(STINGER) +1;
    unique_ptr<uint64_t[]> ptr_result { new uint64_t[N] }; uint64_t* __restrict result = ptr_result.get();

    #pragma omp parallel for
    for(uint64_t v = 0; v < N; v++){
        // timeout check
        if(tcheck->is_timeout()) continue; // do not do any additional processing

        double coeff = 0.0; // by definition, if the vertex has less than two neighbours, its clustering coefficient is zero
        uint64_t degree = stinger_degree_get(STINGER, v);
        if(degree >= 2){
            uint64_t num_triangles = lcc_count_triangles(v);
            uint64_t max_num_edges = degree * (degree -1);
            coeff = static_cast<double>(num_triangles) / max_num_edges;
        }

        result[v] = coeff;
    }

    timer.stop();

    if(tcheck->is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    save(result, N, dump2file);
}

uint64_t StingerDV::lcc_count_triangles(int64_t vertex_id){
    assert(stinger_vtype_get(STINGER, vertex_id) == 0 && "Node marked for deletion");
    uint64_t out = 0; // final number of triangles
    int64_t degree = stinger_degree_get(STINGER, vertex_id);

    // sorted list of neighbours (from both incoming and outgoing edges)
    int64_t* neighbours = (int64_t*) xmalloc(degree * sizeof(int64_t));
    size_t returned_degree { 0 };
    stinger_gather_neighbors(STINGER, vertex_id, &returned_degree, neighbours, /* optional fields weights, timestamps 0 and 1, type */ nullptr, nullptr, nullptr, nullptr, degree);
    assert(returned_degree == static_cast<size_t>(degree) && "Degree mismatch");
    qsort(neighbours, degree, sizeof(int64_t), compare);

    // consider both outgoing and incoming edges. Stinger ensures that a neighbour in a directed
    // graph that has both an incoming and outgoing edge is visited only once
    STINGER_FORALL_EDGES_OF_VTX_BEGIN(STINGER, vertex_id) {
      if (STINGER_EDGE_DEST != vertex_id) { // ignore self loops
        out += lcc_count_intersections(vertex_id, STINGER_EDGE_DEST, neighbours, degree);
      }
    } STINGER_FORALL_EDGES_OF_VTX_END();

    free(neighbours);

    return out;
}

uint64_t StingerDV::lcc_count_intersections (int64_t vertex1, int64_t vertex2, int64_t* vertex1_neighbours, int64_t vertex1_neighbours_sz){
    assert(vertex1_neighbours_sz > 0 && "If this vertex1 doesn't have any neighbours, then vertex2 cannot be a neighbour itself");

    uint64_t out { 0 };

    auto binsearch = [&, vertex1_neighbours, vertex1_neighbours_sz](int64_t vertex3){
        int64_t first = 0;
        int64_t last = vertex1_neighbours_sz -1;
        int64_t middle = (first + last)/2;

        while (first <= last) {
            if (vertex1_neighbours[middle] < vertex3) {
                first = middle + 1;
            } else if (vertex1_neighbours[middle] == vertex3) {
                return true;
            } else {
                last = middle - 1;
            }

            middle = (first + last)/2;
        }

        return false;
    };

    /*
     * According to the Graphalytics spec. v1.0 of this algorithm, for both undirected & directed graphs,
     * we only account the outgoing edges at this step, that is for the connection vertex2 -> vertex3.
     * Whereas for the link vertex1 <-> vertex2 we need to consider both incoming & outgoing edges for
     * directed graphs.
     */
    STINGER_FORALL_OUT_EDGES_OF_VTX_BEGIN(STINGER, vertex2) {
      if (STINGER_EDGE_DEST != vertex1) {
          out += ( binsearch(STINGER_EDGE_DEST) == true );
      }
    } STINGER_FORALL_OUT_EDGES_OF_VTX_END();

    return out;
}


/******************************************************************************
 *                                                                            *
 *  Pagerank                                                                  *
 *                                                                            *
 ******************************************************************************/
// Implementation based on stinger_alg/src/page_rank.c
void StingerDV::pagerank(uint64_t num_iterations, double damping_factor, const char* dump2file){
    auto tcheck = make_unique<utility::TimeoutService>( chrono::seconds{ m_timeout } );
    common::Timer timer; timer.start();

    const size_t num_registered_vertices = stinger_max_active_vertex(STINGER) +1; // logical vertex
    const size_t num_active_vertices = num_vertices();
    auto ptr_rank0 = make_unique<double[]>(num_registered_vertices); double* rank0t = ptr_rank0.get();
    auto ptr_rank1 = make_unique<double[]>(num_registered_vertices); double* rank1t = ptr_rank1.get();

    // init
    #pragma omp parallel for
    for(uint64_t v = 0; v < num_registered_vertices; v++){
        if(stinger_vtype_get(STINGER, v) == 1) continue; // vertex marked for deletion
        rank0t[v] = 1.0 / num_active_vertices;
    }

    for(uint64_t i = 0; i < num_iterations; i++){
        if(tcheck->is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }
        COUT_DEBUG("iteration: " << i);
        double* __restrict rank0 = rank0t; // previous iteration
        double* __restrict rank1 = rank1t; // current iteration

        double pr_constant = 0.0; // dangling sum
        #pragma omp parallel for reduction(+:pr_constant)
        for(uint64_t v = 0; v < num_registered_vertices; v++){
            rank1[v] = 0.0;

            // is this vertex a sink? -> add up to the dangling sum
            if(stinger_outdegree(STINGER, v) == 0){ pr_constant += rank0[v]; }

            // compute the score for the current iteration looking at the incoming edges
            STINGER_FORALL_IN_EDGES_OF_VTX_BEGIN(STINGER, v) {
                uint64_t outdegree = stinger_outdegree (STINGER, STINGER_EDGE_DEST);
                if(outdegree > 0){ rank1[v] += rank0[STINGER_EDGE_DEST] / outdegree; }
            } STINGER_FORALL_IN_EDGES_OF_VTX_END();
        }

        #pragma omp parallel for
        for(uint64_t v = 0; v < num_registered_vertices; v++){
            rank1[v] = (1.0 - damping_factor) / num_active_vertices + (rank1[v] + (pr_constant / num_active_vertices)) * damping_factor;
        }

        std::swap(rank0t, rank1t); // next iteration
    }

    if(tcheck->is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // store the final results (if required)
    save(rank0t, num_registered_vertices, dump2file);
};

/******************************************************************************
 *                                                                            *
 *  Shortest paths                                                            *
 *                                                                            *
 *****************************************************************************/
// copy & paste from stinger_alg/src/shortest_paths.cpp
namespace {
typedef struct {
    int64_t vertex;
    double  cost;
} weighted_vertex_t;

//this is a simple comparison function for the queue that sorts based on the cost to reach a vertex
static bool comp(weighted_vertex_t a, weighted_vertex_t b) {
    return a.cost > b.cost;
}

typedef bool(*CompareType)(weighted_vertex_t a, weighted_vertex_t b);

static std::vector<double> a_star(stinger_t * S, int64_t NV, int64_t source_vertex, int64_t dest_vertex, bool ignore_weights, chrono::seconds timeout){
    auto tcheck = make_unique<utility::TimeoutService>( timeout );
    common::Timer timer; timer.start();
    std::vector<double> cost_so_far(NV);

    for (int64_t v = 0; v < NV; v++){
        cost_so_far[v] = std::numeric_limits<double>::max(); // initialize all the values to + infinity
    }

    //represents the vertices to look at in the graph
    std::priority_queue<weighted_vertex_t, std::vector<weighted_vertex_t>, CompareType> frontier (comp);
    //add the initial vertex to the priority queue
    weighted_vertex_t source;
    source.vertex = source_vertex;
    source.cost = 0;
    cost_so_far[source_vertex] = 0;
    frontier.push(source);

    while (!frontier.empty()) {
        if(tcheck->is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

        weighted_vertex_t current = frontier.top();
        frontier.pop();

        //we found our goal!
        if (current.vertex == dest_vertex) break;

        STINGER_FORALL_OUT_EDGES_OF_VTX_BEGIN(S, current.vertex) {
            //for all the neighbors of the current vertex
            double new_cost;
            if (ignore_weights){
                new_cost = cost_so_far[current.vertex] + 1; // disregard the weights, instead count paths
            } else{
                new_cost = cost_so_far[current.vertex] + INT2DBL(STINGER_EDGE_WEIGHT);
            }

            if (new_cost < cost_so_far[STINGER_EDGE_DEST] || cost_so_far[STINGER_EDGE_DEST] == std::numeric_limits<double>::max()) {
                cost_so_far[STINGER_EDGE_DEST] = new_cost;
                weighted_vertex_t next;
                next.vertex = STINGER_EDGE_DEST;
                next.cost = new_cost;
                frontier.push(next);
            }
        }
        STINGER_FORALL_OUT_EDGES_OF_VTX_END();
    }

    return cost_so_far;
}

// dump the content to the given file
static void save_shortest_paths(const vector<double>& result, bool weighted, const char* dump2file){
    if(dump2file == nullptr) return; // nop
    COUT_DEBUG("save the results to: " << dump2file)
    string label_infinity = weighted ? string("infinity") : to_string(numeric_limits<uint64_t>::max());

    fstream handle(dump2file, ios_base::out);
    if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

    for(uint64_t vertex_id = 0; vertex_id < result.size(); vertex_id++){
        handle << vertex_id << " ";
        if(result[vertex_id] == numeric_limits<double>::max()){
            handle << label_infinity;
        } else {
            handle << result[vertex_id];
        }
        handle << "\n";
    }

    handle.close();
}
} // anonymous namespace

void StingerDV::bfs(uint64_t source_vertex_id, const char* dump2file){
    auto result = a_star(STINGER, stinger_max_active_vertex(STINGER) +1, source_vertex_id, /* all vertices */ -1, /* ignore weights */ true, chrono::seconds( m_timeout ));
    save_shortest_paths(result, /* weighted ? */ false, dump2file);
}

void StingerDV::sssp(uint64_t source_vertex_id, const char* dump2file){
    auto result = a_star(STINGER, stinger_max_active_vertex(STINGER) +1, source_vertex_id, /* all vertices */ -1, /* ignore weights */ false, chrono::seconds( m_timeout ));
    save_shortest_paths(result, /* weighted ? */ true, dump2file);
}

/******************************************************************************
 *                                                                            *
 *  Dump                                                                      *
 *                                                                            *
 *****************************************************************************/
void StingerDV::dump_ostream(ostream& out) const {
    struct stinger* graph = STINGER;

    out << "[STINGER] Vertices: " << num_vertices() << ", edges: " << num_edges() << ", directed: " << std::boolalpha << is_directed() << ", size: " << stinger_graph_size(graph) << " bytes" << "\n";
    for(int64_t vertex_id = 0, vertex_max = stinger_max_active_vertex(graph); vertex_id <= vertex_max; vertex_id++){

        out << "[" << vertex_id << "] degree in/out/total: " << stinger_indegree_get(graph, vertex_id) << "/" << stinger_outdegree_get(graph, vertex_id) << "/" << stinger_degree_get(graph, vertex_id);
        if(stinger_outdegree_get(graph, vertex_id) == 0){
            out << "\n";
        } else {
            out << ", outgoing edges: \n";
            STINGER_FORALL_OUT_EDGES_OF_VTX_BEGIN(graph, vertex_id) {
                out << "  " << STINGER_EDGE_DEST << ", weight: " << INT2DBL(STINGER_EDGE_WEIGHT) << "\n";
            } STINGER_FORALL_OUT_EDGES_OF_VTX_END();
        }
    }
    std::flush(out);
}

template<typename T>
void StingerDV::save(const T& data, uint64_t sz, const char* dump2file){
    if(dump2file == nullptr) return; // nop
    COUT_DEBUG("save the results to: " << dump2file)

    fstream handle(dump2file, ios_base::out);
    if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

    for(uint64_t v = 0; v < sz; v++){
        handle << v << " " << data[v] << "\n";
    }

    handle.close();
}


} // namespace
