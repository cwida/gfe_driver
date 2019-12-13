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

#include "graphone.hpp"
#include "internal.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <omp.h>
#include <unordered_map>
#include <unordered_set>

#include "common/system.hpp"
#include "common/timer.hpp"
#include "third-party/libcuckoo/cuckoohash_map.hh"
#include "utility/timeout_service.hpp"

using namespace common;
using namespace gfe::utility;
using namespace std;

/*****************************************************************************
 *                                                                           *
 *  Globals                                                                  *
 *                                                                           *
 *****************************************************************************/
// these globals are required by the current implementation of the GraphOne library...
graph* g = nullptr;
int THD_COUNT = 0;

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG
namespace gfe { extern mutex _log_mutex [[maybe_unused]]; }
#define COUT_DEBUG_FORCE(msg) { std::scoped_lock<std::mutex> lock{::gfe::_log_mutex}; std::cout << "[GraphOne::" << __FUNCTION__ << "] [Thread #" << common::concurrency::get_thread_id() << "] " << msg << std::endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif

/*****************************************************************************
 *                                                                           *
 *  Init                                                                     *
 *                                                                           *
 *****************************************************************************/

namespace gfe::library {

GraphOne::GraphOne(bool is_graph_directed, bool use_vertex2id_mapping, bool blind_writes, bool ignore_build, uint64_t max_num_vertices) :
        m_is_directed(is_graph_directed), m_translate_vertex_ids(use_vertex2id_mapping), m_blind_writes(blind_writes), m_ignore_build(ignore_build) {
    if(g != nullptr) ERROR("An instance of GraphOne has already been created");

    // The graph instance is a global ...
    ::g = new class graph();

    // THD_COUNT is another global
    ::THD_COUNT = (omp_get_max_threads() -1); // as defined by the experiments in the GraphOne suite

    g->cf_info = new cfinfo_t*[2]; // number of containers ("column families") = 2: 1. metadata (vertex dictionary), 2. edges of the graph
    g->p_info = nullptr; // property names, unused

    // metadata
    auto metadata = new typekv_t(); // vertex types or classes
    g->add_columnfamily(metadata); // cf_info[0] <- typekv_t

    // the actual weighted graph
    pgraph_t<lite_edge_t>* weighted_graph { nullptr };
    if(!m_is_directed){
        weighted_graph = new p_ugraph_t(); // undirected
    } else {
        weighted_graph = new p_dgraph_t(); // directed
    }
    weighted_graph->flag1 = weighted_graph->flag2 =1; // `1' is the vertex type
    g->add_columnfamily(weighted_graph); // register the graph in the database instance

    metadata->manual_setup(/* max number of vertices */ max_num_vertices, /* statically created the vertices ? */ false);
    g->prep_graph_baseline(); // I think this finalises the schema
    g->create_threads(/* archiving */ true, /* logging */ false); // background threads, to create the snapshots and logging to disk

    // edge locks, to ensure consistency when performing an update
    if(use_vertex2id_mapping || !blind_writes){
        m_num_edge_locks = 64; // arbitrary value
        static_assert(sizeof(PaddedLock) == 64, "Expected to match the size of a cache line");
        m_edge_locks = new PaddedLock[ m_num_edge_locks ]();
    }
}

GraphOne::~GraphOne(){
    delete[] m_edge_locks; m_edge_locks = nullptr; m_num_edge_locks = 0;

    delete ::g; ::g = nullptr;
}

/*****************************************************************************
 *                                                                           *
 *  Properties                                                               *
 *                                                                           *
 *****************************************************************************/

bool GraphOne::is_directed() const {
    return m_is_directed;
}

bool GraphOne::is_undirected() const {
    return !m_is_directed;
}

uint64_t GraphOne::num_vertices() const {
    return m_num_vertices;
}

uint64_t GraphOne::num_edges() const {
    return m_num_edges;
}

bool GraphOne::has_blind_writes() const {
    return m_blind_writes;
}

bool GraphOne::has_vertex(uint64_t vertex_id) const {
    if(m_translate_vertex_ids){
        string str_vertex_id = to_string(vertex_id);
        sid_t logic_vertex_id = g->get_typekv()->get_sid(str_vertex_id.c_str());
        return logic_vertex_id != INVALID_SID;
    } else {
        return vertex_id < m_num_vertices;
    }
}

uint64_t GraphOne::vtx_ext2int(uint64_t external_vertex_id) const {
    if(m_translate_vertex_ids){
        string str_source_vertex_id = to_string(external_vertex_id);
        sid_t sid = g->get_sid(str_source_vertex_id.c_str());
        if(sid == INVALID_SID) ERROR("[1] Invalid vertex ID: " << external_vertex_id);
        return static_cast<uint64_t>(sid);
    } else {
        if(external_vertex_id >= m_num_vertices) ERROR("[2] Invalid vertex ID: " << external_vertex_id);
        return external_vertex_id;
    }
}

bool GraphOne::find_edge(uint64_t v0, uint64_t v1, bool* out_record_found, double* out_weight) const {
    auto* view = create_static_view(get_graphone_graph(), /* default mask */ 0); // default = adjacency list of the latest snapshot + non archived edges
    weight_edge_t* edges { nullptr };
    int64_t num_edges = (int64_t) view->get_nonarchived_edges(edges);
    bool found = false; bool exists = false;

    // first check the non archived edges, that is, the write store
    int64_t i = num_edges -1;
    while(i >= 0 && !found){
        uint64_t src = TO_VID(get_src(edges[i]));
        uint64_t dst = get_dst(edges[i]);
        bool is_deletion = IS_DEL(get_src(edges[i]));

        if(src == v0 && dst == v1){
            found = true; exists = !is_deletion;
            if(out_weight != nullptr) { *out_weight = edges[i].dst_id.second.value_double; }
        } else {
            i--;
        }
    }

    // search among the archived edges, that is, the read store
    if(!found){
        unique_ptr<lite_edge_t[]> ptr_neighbours { nullptr };
        weight_sid_t* neighbours = nullptr;
        uint64_t neighbours_sz = 0;

        // we want to iterate on the adjacency list with less neighbours, whether this is
        // represented by the outgoing edges of v0 or the incoming edges of v1
        uint64_t degree_v0 = view->get_degree_out(v0);
        uint64_t degree_v1 = view->get_degree_in(v1);
        uint64_t vertex_to_match = -1;

        // access the adjacency list with the minor degree
        if(degree_v0 <= degree_v1){
            neighbours_sz = degree_v0;
            ptr_neighbours.reset( new lite_edge_t[neighbours_sz] );
            neighbours = ptr_neighbours.get();
            // it turns out that degree_v0 and degree_v1 are only upper bounds, we get the actual degree only when we acquire the edges
            view->get_nebrs_out(v0, neighbours);
            vertex_to_match = v1;
        } else {
            neighbours_sz = degree_v1;
            ptr_neighbours.reset( new lite_edge_t[neighbours_sz] );
            neighbours = ptr_neighbours.get();
            view->get_nebrs_in(v1, neighbours);
            vertex_to_match = v0;
        }

//        COUT_DEBUG("find_edge " << v0 << " -> " << v1 << ", degree_v0: " << degree_v0 << ", degree_v1: " << degree_v1 << ", neighbours_sz: " << neighbours_sz);

        // search the target match, either v0 or v1, depending on which adjacency list has been selected
        i = static_cast<int64_t>(neighbours_sz) - 1;
        while(i >= 0 && !found){
            uint64_t vertex_id = TO_VID(get_sid(neighbours[i]));
            bool is_deletion = IS_DEL(get_sid(neighbours[i]));
            assert(is_deletion == false && "I think in the current implementation of GraphOne, the adj list only returns the non deleted vertices");
            if(vertex_id == vertex_to_match){
                found = true; exists = !is_deletion;
                if(out_weight != nullptr) { *out_weight = neighbours[i].second.value_double; }
            } else {
                i--;
            }
        }
    }

    delete_static_view(view);

    if(out_record_found != nullptr) *out_record_found = found;
    return exists;
}

double GraphOne::get_weight(uint64_t source, uint64_t destination) const {
    COUT_DEBUG("source: " << source << ", destination: " << destination);

    auto nan = numeric_limits<double>::quiet_NaN();

    sid_t v0, v1;
    if(m_translate_vertex_ids){
        string str_source = to_string(source);
        string str_destination = to_string(destination);

        auto labels = g->get_typekv();
        v0 = labels->get_sid(str_source.c_str());
        if(v0 == INVALID_SID) return nan;
        v1 = labels->get_sid(str_destination.c_str());
        if(v1 == INVALID_SID) return nan;
    } else {
        v0 = (sid_t) source;
        v1 = (sid_t) destination;
    }

    if(!m_is_directed && v0 > v1){ std::swap(v0, v1); }

    return get_weight_impl(v0, v1);
}

double GraphOne::get_weight_impl(uint64_t v0, uint64_t v1) const {
    double weight = 0.0;
    bool exists = find_edge(v0, v1, /* record found ? */ nullptr, &weight);
    return exists ? weight : numeric_limits<double>::quiet_NaN();
}

void GraphOne::set_timeout(uint64_t seconds){
    m_timeout = chrono::seconds{ seconds };
}

/*****************************************************************************
 *                                                                           *
 *  Updates                                                                  *
 *                                                                           *
 *****************************************************************************/
bool GraphOne::add_vertex(uint64_t vertex_id) {
    COUT_DEBUG("Vertex: " << vertex_id);

    if(m_translate_vertex_ids){ // regular graph, translate the vertex into the logical id into the adjacency list
        string str_vertex_id = to_string(vertex_id);

//        auto& mutex = m_edge_locks[vertex_id % m_num_edge_locks].m_lock;
        // The internal dictionary provided by GraphOne is not thread safe
        scoped_lock<SpinLock> lock(m_mutex_vtx);

        sid_t id = g->type_update(str_vertex_id);
        bool vertex_exists = ( id == INVALID_SID );

        if(!vertex_exists){
            COUT_DEBUG("vertex_id " << vertex_id << ", new id: " << id);
            m_num_vertices++;
            return true;
        } else { // the vertex already exists
            return false;
        }

    } else { // dense graph, without the translation map
        uint64_t num_vertices = m_num_vertices;
        uint64_t new_value = std::max<uint64_t>(num_vertices, vertex_id +1);

        while(!m_num_vertices.compare_exchange_weak(/* by ref, out */ num_vertices, /* by value */ new_value,
                /* memory order in case of success */ std::memory_order_release,
                /* memory order in case of failure */ std::memory_order_relaxed)){
            new_value = std::max<uint64_t>(num_vertices, vertex_id +1);
        }

        return vertex_id +1 == new_value;
    }
}

bool GraphOne::remove_vertex(uint64_t vertex_id){
    COUT_DEBUG("Vertex: " << vertex_id);

//    if(m_translate_vertex_ids){ // regular graph, translate the vertex into the logical id into the adjacency list
//        string str_vertex_id = to_string(vertex_id);
//        auto labels = g->get_typekv();
//
//        scoped_lock<SpinLock> lock(m_mutex);
//        sid_t logic_vertex_id = labels->str2vid.find(str_vertex_id);
//        if(logic_vertex_id == INVALID_SID){
//            return false; // the mapping does not exist
//        } else {
//            labels->str2vid.erase(str_vertex_id);
//            g->vert_count--;
//            return true;
//        }
//
//    } else { // dense graph, vertices cannot be removed
//        assert(m_num_vertices > 0);
//        m_num_vertices--;
//
//        // ...
//        return false;
//    }

    // vertices cannot be removed atm...
    return false;
}

bool GraphOne::add_edge(gfe::graph::WeightedEdge e){
    COUT_DEBUG("Edge: " << e);

    sid_t v0 {0}, v1 {0};

    if(m_translate_vertex_ids){
        auto labels = g->get_typekv();

        string str_source = to_string(e.source());
        v0 = labels->get_sid(str_source.c_str());
        if(v0 == INVALID_SID) return false; // the source vertex does not exist

        string str_destination = to_string(e.destination());
        v1 = labels->get_sid(str_destination.c_str());
        if(v1 == INVALID_SID) return false; // the destination vertex does not exist
    } else {
        v0 = (sid_t) e.source();
        v1 = (sid_t) e.destination();
    }

    // not strictly necessary, but it eases the impl of #get_weight
    if(!m_is_directed && v0 > v1) std::swap(v0, v1);

    COUT_DEBUG("v0: " << v0 << ", v1: " << v1);

    if(has_blind_writes()){
        do_update(/* is insert ? */ true, v0, v1, e.weight());
        return true;
    } else {
        SpinLock& mutex = m_edge_locks[(v0 + v1) % m_num_edge_locks].m_lock;
        scoped_lock<SpinLock> xlock(mutex);

        bool exists = find_edge(v0, v1);
        if(!exists){
            do_update(/* is insert ? */ true, v0, v1, e.weight());
            return true;
        } else {
            return false;
        }
    }
}

bool GraphOne::remove_edge(gfe::graph::Edge e){
    COUT_DEBUG("Edge: " << e);

    sid_t v0 {0}, v1 {0};

    if(m_translate_vertex_ids){
        auto labels = g->get_typekv();

        string str_source = to_string(e.source());
        v0 = labels->get_sid(str_source.c_str());
        if(v0 == INVALID_SID) return false; // the source vertex does not exist

        string str_destination = to_string(e.destination());
        v1 = labels->get_sid(str_destination.c_str());
        if(v1 == INVALID_SID) return false; // the destination vertex does not exist
    } else {
        v0 = (sid_t) e.source();
        v1 = (sid_t) e.destination();
    }

    // not strictly necessary, but it eases the impl of #get_weight
    if(!m_is_directed && v0 > v1) std::swap(v0, v1);

    if(has_blind_writes()){
        do_update(/* is insert ? */ false, v0, v1);
        return true;
    } else {
        SpinLock& mutex = m_edge_locks[(v0 + v1) % m_num_edge_locks].m_lock;
        scoped_lock<SpinLock> xlock(mutex);

        bool exists = find_edge(v0, v1);
        if(exists){
            do_update(/* is insert ? */ false, v0, v1);
            return true;
        } else {
            return false;
        }
    }
}

void GraphOne::do_update(bool is_insert, uint64_t v0, uint64_t v1, double weight){
    assert((is_directed() || (v0 < v1)) && "At least in raw insertions in the archive buffer, edges need to be inserted with src < dst");
    COUT_DEBUG("is_insert: " << is_insert << ", v0: " << v0 << ", v1: " << v1);

    weight_edge_t edge;
    if(is_insert){
        set_src(edge, v0); // edge source
    } else {
        set_src(edge, DEL_SID(v0)); // edge source + mask to ask for deletion
    }
    set_dst(edge, v1); // destination
    edge.dst_id.second.value_double = weight;

    // Perform the update
    status_t rc = get_graphone_graph()->batch_edge(edge);
    if(rc == eEndBatch) m_num_levels++;

    // update the global counter on the number of edges present
    if(is_insert){
        m_num_edges++; // atomic
    } else {

        assert(m_num_edges > 0);
        m_num_edges--; // atomic
    }
}

void GraphOne::build(){
    COUT_DEBUG("Build");
    if(m_ignore_build) return; // nop
    g->waitfor_archive();
    m_num_levels++;
}

uint64_t GraphOne::num_levels() const {
    return m_num_levels;
}

/*****************************************************************************
 *                                                                           *
 *  Dump                                                                     *
 *                                                                           *
 *****************************************************************************/

void GraphOne::dump_ostream(std::ostream& out) const {
    uint64_t N = num_vertices();
    auto view = create_static_view(get_graphone_graph(), STALE_MASK);

    out << "[GraphOne] directed graph: " << boolalpha << is_directed() << ", use vertex "
            "dictionary: " << m_translate_vertex_ids << ", blind writes: " << has_blind_writes() << ", num vertices: " << N << ", vertex array capacity: " << view->get_vcount() << ", num edges: " << num_edges() << endl;

    // only consider the archived edges
    lite_edge_t* neighbours = nullptr;
    uint64_t neighbours_sz = 0;

    for(sid_t vertex_id = 0; vertex_id < N; vertex_id ++){
        if(m_translate_vertex_ids){
            string vertex_name = g->get_typekv()->get_vertex_name(vertex_id);
            if(vertex_name.empty()) continue;
            out << "[vertex_id: " << TO_VID(vertex_id) << ", name: " << vertex_name << "]\n";
        } else {
            out << "[vertex_id: " << TO_VID(vertex_id) << "]\n";
        }

        uint64_t degree_out = view->get_degree_out(vertex_id); // it used to be an upper bound, the actual degree may be changed by view->get_nebrs_out(...)
        if(degree_out > neighbours_sz){
            neighbours_sz = degree_out;
            neighbours = (lite_edge_t*) realloc(neighbours, sizeof(neighbours[0]) * degree_out);
            if(neighbours == nullptr) throw std::bad_alloc{};
        }
        degree_out = view->get_nebrs_out(vertex_id, neighbours);
        out << "  " << degree_out << " outgoing edges: ";

        for(uint64_t edge_id = 0; edge_id < degree_out; edge_id++){
            if(edge_id > 0) cout << ", ";
            cout << get_sid(neighbours[edge_id]) << " [w=" << neighbours[edge_id].second.value_double << "]";
        }
        out << "\n";

        if(is_directed()){
            uint64_t degree_in = view->get_degree_in(vertex_id); // it used to be an upper bound, the actual degree may be changed by view->get_nebrs_in(...)
            if(degree_in > neighbours_sz){
                neighbours_sz = degree_out;
                neighbours = (lite_edge_t*) realloc(neighbours, sizeof(neighbours[0]) * degree_in);
                if(neighbours == nullptr) throw std::bad_alloc{};
            }
            degree_in = view->get_nebrs_in(vertex_id, neighbours);
            out << "  " << degree_in << " incoming edges: ";
            for(uint64_t edge_id = 0; edge_id < degree_in; edge_id++){
                if(edge_id > 0) cout << ", ";
                out << get_sid(neighbours[edge_id]) << " [w=" << neighbours[edge_id].second.value_double << "]";
            }

            out << "\n";
        }
    }

    free(neighbours); neighbours = nullptr; neighbours_sz = 0;

    weight_edge_t* edges2 { nullptr };
    index_t edges2_sz = view->get_nonarchived_edges(edges2);
    if(edges2_sz == 0){
        out << "There are no non-archived edges" << endl;
    } else {
        out << "There are " << edges2_sz << " non-archived edges: ";

        for (index_t i = 0; i < edges2_sz; ++i) {
            if(i > 0) cout << ", ";
            sid_t v0 = get_src(edges2[i]);
            sid_t v1 = get_dst(edges2[i]);
            double weight = edges2[i].dst_id.second.value_double;
            cout << "< " << v0 << ", " << v1 << " [" << weight << "]" << ">";
        }
    }

    delete_static_view(view);
}

/*****************************************************************************
 *                                                                           *
 *  BFS                                                                      *
 *                                                                           *
 *****************************************************************************/
// Implementation based on:
// > file graphone/test/plaingraph_test.cpp
// >>> 03/12/2019, the function was moved to analytics/mem_iterative_analytics.h after the latest upstream bump
// > function: template<class T> void mem_bfs(gview_t<T>* snaph, uint8_t* status, sid_t root)

static unique_ptr<uint32_t[]> graphone_execute_bfs(uint64_t v_count, uint64_t root, utility::TimeoutService& timeout){
    // Initialisation
    auto* view = create_static_view(get_graphone_graph(), SIMPLE_MASK); // global
    int level = 1;
    int top_down = 1;
    sid_t frontier = 0;
//    sid_t v_count = num_vertices(); // snaph->get_vcount(); // parameter
    unique_ptr<uint32_t[]> ptr_status { new uint32_t[v_count]() }; // avoid memory leaks
    uint32_t* status = ptr_status.get();
    status[root] = level; // 0 => unreachable vertices, 1 = root, >1 = all the other reachable vertices

    // BFS Algorithm
    do {
        frontier = 0;

        #pragma omp parallel reduction(+:frontier)
        {
            // the array with all the edges of a given vertex
            lite_edge_t* neighbours = nullptr;
            uint64_t neighbours_sz = 0;

            if (top_down) {

                #pragma omp for nowait
                for (vid_t v = 0; v < v_count; v++) {
                    if (status[v] != level) continue;
                    uint64_t degree = view->get_degree_out(v);
                    if(degree == 0) continue;

                    // The original code from
                    // makes use of view->get_nebrs_archived_out(), but t only retrieves the change sets, without
                    // merging the deletions. Hence, that code assumes that only insertions occurred.

                    if(degree > neighbours_sz){
                        neighbours_sz = degree;
                        neighbours = (lite_edge_t*) realloc(neighbours, sizeof(neighbours[0]) * neighbours_sz);
                        if(neighbours == nullptr) throw std::bad_alloc{};
                    }
                    degree = view->get_nebrs_out(v, neighbours); // reassign degree as it may have been decreased by #get_nebrs_out

                    for(uint64_t i = 0; i < degree; i++){
                        uint64_t dst = get_sid(neighbours[i]);
                        if(status[dst] == 0){ // this node has not been visited yet
                            status[dst] = level +1;
                            ++frontier;
                        }
                    }
                }

            } else { //bottom up

                #pragma omp for nowait
                for (vid_t v = 0; v < v_count; v++) {
                    if (status[v] != 0) continue; // this node has already been visited
                    uint64_t degree = view->get_degree_in(v);
                    if (0 == degree) continue;

                    if(degree > neighbours_sz){
                        neighbours_sz = degree;
                        neighbours = (lite_edge_t*) realloc(neighbours, sizeof(neighbours[0]) * neighbours_sz);
                        if(neighbours == nullptr) throw std::bad_alloc{};
                    }

                    degree = view->get_nebrs_in(v, neighbours); // reassign degree as it may have been decreased by #get_nebrs_in

                    for(uint64_t i = 0; i < degree; i++){
                        uint64_t src = get_sid(neighbours[i]);
                        if(status[src] == level){
                            status[v] = level +1;
                            ++frontier;
                            break;
                        }
                    }
                }
            } // end if

            free(neighbours); neighbours = nullptr; neighbours_sz = 0;
        }

        //Point is to simulate bottom up bfs, and measure the trade-off
        if ((frontier >= 0.002*v_count) && ( 0 == view->is_unidir())){ // || level == 2)
            top_down = false;
        } else {
            top_down = true;
        }
        ++level;

    } while (frontier && !timeout.is_timeout());

    delete_static_view(view);

    return ptr_status;
}

void GraphOne::bfs(uint64_t source_vertex_id, const char* dump2file) {
    // Init
    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();
    uint64_t root = vtx_ext2int(source_vertex_id); // from a graphalytics vertex ID to the dense vertex ID used internally
    uint64_t N = num_vertices();

    // Run the BFS algorithm
    auto ptr_distances = graphone_execute_bfs(N, root, timeout);
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);  }

    // Translate the vertex IDs and dump to file if required
    uint32_t* distances = ptr_distances.get();
    if(m_translate_vertex_ids) { // translate from GraphOne vertex ids to external vertex ids
        cuckoohash_map</* external id */ string, /* distance */ int> external_ids;
#pragma omp parallel for
        for(sid_t internal_id = 0; internal_id < N; internal_id++){
            // 1.what is the real node ID, in the external domain (e.g. user id)
            string vertex_name = g->get_typekv()->get_vertex_name(internal_id);
            if(vertex_name.empty()) continue; // this means that a mapping does not exist. It should never occur as atm we don't support vertex deletions

            // 2. retrieve the distance / weight
            int distance = distances[internal_id];

            // 3. make the association vertex name - distance
            external_ids.insert(vertex_name, distance);
        }

        if(timeout.is_timeout()){
            RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);
        }

        // store the results in the given file
        if(dump2file != nullptr){
            COUT_DEBUG("save the results to: " << dump2file);
            fstream handle(dump2file, ios_base::out);
            if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

            auto hashtable = external_ids.lock_table();

            for(const auto& keyvaluepair : hashtable){
                handle << keyvaluepair.first << " ";
                auto distance = keyvaluepair.second;

                // distance = 0 => node never visited
                // distance = 1 => root
                // distance > 1 => other nodes

                if(distance > 0){
                    handle << distance -1;
                } else {
                    handle << std::numeric_limits<int64_t>::max(); // it should have been -1, but ok
                }

                handle << "\n";
            }

            handle.close();
        }

    } else { // without the vertex dictionary

        // store the results in the given file
        if(dump2file != nullptr){
            COUT_DEBUG("save the results to: " << dump2file);
            fstream handle(dump2file, ios_base::out);
            if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");


            for(sid_t internal_id = 0; internal_id < N; internal_id++){
                handle << internal_id << " ";
                auto distance = distances[internal_id];

                // distance = 0 => node never visited
                // distance = 1 => root
                // distance > 1 => other nodes

                if(distance > 0){
                    handle << distance -1;
                } else {
                    handle << std::numeric_limits<int64_t>::max(); // it should have been -1, but ok
                }

                handle << "\n";
            }

            handle.close();
        }
    }
}

/*****************************************************************************
 *                                                                           *
 *  PageRank                                                                 *
 *                                                                           *
 *****************************************************************************/
// GraphOne ships with its own implementation of PageRank in graphone/analytics/mem_iterative_analytics.h,
// function template<class T> void mem_pagerank(gview_t<T>* snaph, int iteration_count). Still
// there are a few issues w.r.t. to the Graphalytics specification:
// 1. Similarly to the BFS snippet, it operates on change sets and does not handle deletions
// 2. The dangling sum is computed once at the start, differently from what expected by the Graphalytics spec
// 3. It operates on floats rather than doubles
// 4. It uses a different formula to compute the score at the end
// The algorithm is too different, the implementation below is adapted from llama_pagerank.cpp:
//#define DEBUG_PAGERANK
#if defined(DEBUG_PAGERANK)
#define COUT_DEBUG_PAGERANK(msg) COUT_DEBUG(msg)
#else
#define COUT_DEBUG_PAGERANK(msg)
#endif

static unique_ptr<double[]> graphone_execute_pagerank(uint64_t num_vertices, uint64_t num_iterations, double d, utility::TimeoutService& timer){
    COUT_DEBUG_PAGERANK("num_vertices: " << num_vertices << ", num_iterations: " << num_iterations << ", damping factor: " << d);

    auto* view = create_static_view(get_graphone_graph(), SIMPLE_MASK); // global

    unique_ptr<double[]> ptr_G_pg_rank{ new double[num_vertices]() }; // avoid memory leaks
    unique_ptr<double[]> ptr_G_pg_rank_nxt{ new double[num_vertices]() };
    double* G_pg_rank = ptr_G_pg_rank.get();
    double* G_pg_rank_nxt = ptr_G_pg_rank_nxt.get();
    double N = num_vertices;

    // init
    #pragma omp parallel for
    for(uint64_t v = 0; v < num_vertices; v++){
        G_pg_rank[v] = 1.0 / N;
        G_pg_rank_nxt[v] = 0.0;
    }

    // iterations
    uint64_t current_iteration = 0;
    while(current_iteration < num_iterations && !timer.is_timeout()){

        double dsum = 0.0; // dangling sum
        #pragma omp parallel reduction(+:dsum)
        {
            lite_edge_t* neighbours = nullptr;
            uint64_t neighbours_sz = 0;

            #pragma omp for
            for(uint64_t v = 0; v < num_vertices; v++){
                uint64_t degree_in = view->get_degree_in(v);
                uint64_t degree_out = view->get_degree_out(v);

                if(degree_out == 0){ // this is a sink
                    dsum += G_pg_rank[v];
                }

                if(degree_in > neighbours_sz){
                    neighbours_sz = degree_in;
                    neighbours = (lite_edge_t*) realloc(neighbours, sizeof(neighbours[0]) * neighbours_sz);
                    if(neighbours == nullptr) throw std::bad_alloc{};
                }
                view->get_nebrs_in(v, neighbours);

                double score = 0;
                for(uint64_t i = 0; i < degree_in; i++){
                    uint64_t u = get_sid(neighbours[i]);
                    score += G_pg_rank[u] / view->get_degree_out(u);
                }

                G_pg_rank_nxt[v] = score;
            }

            free(neighbours); neighbours = nullptr; neighbours_sz = 0;
        }


        // rank <- rank_next; rank_next <- 0;
        #pragma omp for schedule(dynamic, 4096)
        for (uint64_t v = 0; v < num_vertices; v++) {
            G_pg_rank[v] = (1-d) / N + d * G_pg_rank_nxt[v] + /* dangling sum */ d * dsum / N;
            G_pg_rank_nxt[v] = 0.0;
        }

        // next iteration
        current_iteration++;
    }

    delete_static_view(view);

    return ptr_G_pg_rank;
}

void GraphOne::pagerank(uint64_t num_iterations, double damping_factor, const char* dump2file) {
    // Init
    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();
    uint64_t N = num_vertices();

    // Run the PageRank algorithm
    unique_ptr<double[]> ptr_rank = graphone_execute_pagerank(N, num_iterations, damping_factor, timeout);
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);  }

    // Translate the vertex IDs and dump to file if required
    double* rank = ptr_rank.get();
    if(m_translate_vertex_ids) { // translate from GraphOne vertex ids to external vertex ids
        cuckoohash_map</* external id */ string, /* rank */ double> external_ids;
        #pragma omp parallel for
        for(sid_t internal_id = 0; internal_id < N; internal_id++){
            // 1.what is the real node ID, in the external domain (e.g. user id)
            string vertex_name = g->get_typekv()->get_vertex_name(internal_id);
            if(vertex_name.empty()) continue; // this means that a mapping does not exist. It should never occur as atm we don't support vertex deletions

            // 2. retrieve the distance / weight
            double score = rank[internal_id];

            // 3. make the association vertex name - score
            external_ids.insert(vertex_name, score);
        }

        if(timeout.is_timeout()){
            RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);
        }

        // store the results in the given file
        if(dump2file != nullptr){
            COUT_DEBUG("save the results to: " << dump2file);
            fstream handle(dump2file, ios_base::out);
            if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

            auto hashtable = external_ids.lock_table();

            for(const auto& keyvaluepair : hashtable){
                handle << keyvaluepair.first << " " << keyvaluepair.second << "\n";
            }

            handle.close();
        }

    } else { // without the vertex dictionary

        // store the results in the given file
        if(dump2file != nullptr){
            COUT_DEBUG("save the results to: " << dump2file);
            fstream handle(dump2file, ios_base::out);
            if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

            for(sid_t internal_id = 0; internal_id < N; internal_id++){
                handle << internal_id << " " << rank[internal_id] << "\n";
            }

            handle.close();
        }
    }
}

/*****************************************************************************
 *                                                                           *
 *  WCC                                                                      *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG_WCC
#if defined(DEBUG_WCC)
#define COUT_DEBUG_WCC(msg) COUT_DEBUG(msg)
#else
#define COUT_DEBUG_WCC(msg)
#endif

// Implementation similar to the specification
static unique_ptr<uint32_t[]> graphone_execute_wcc(uint64_t num_vertices, bool is_directed, utility::TimeoutService& timeout){
    // Initialisation
    auto* view = create_static_view(get_graphone_graph(), SIMPLE_MASK); // global
    unique_ptr<uint32_t[]> ptr_components0 { new uint32_t[num_vertices] }; // avoid memory leaks
    unique_ptr<uint32_t[]> ptr_components1 { new uint32_t[num_vertices] };
    uint32_t* components0 = ptr_components0.get(); // current iteration
    uint32_t* components1 = ptr_components1.get(); // next iteration

    #pragma omp parallel for
    for(uint64_t v = 0; v < num_vertices; v++){
        components0[v] = v;
    }

    bool converged = false; // shared var
    do {

        converged = true;
        #pragma omp parallel shared(converged)
        {
            lite_edge_t* neighbours = nullptr;
            uint64_t neighbours_sz = 0;

            #pragma omp for
            for(uint64_t v = 0; v < num_vertices; v++){
                uint64_t degree_out = view->get_degree_out(v);
                uint64_t value_next = components0[v]; // init with the old value

                // outgoing edges
                if(degree_out > neighbours_sz){
                    COUT_DEBUG_WCC("realloc the neighbours (out): " << neighbours_sz << " -> " << degree_out);
                    neighbours_sz = degree_out;
                    neighbours = (lite_edge_t*) realloc(neighbours, sizeof(neighbours[0]) * neighbours_sz);
                    if(neighbours == nullptr) throw std::bad_alloc{};
                }
                view->get_nebrs_out(v, neighbours);

                for(uint64_t i = 0; i < degree_out; i++){
                    uint64_t u = get_sid(neighbours[i]);
                    if(components0[u] < value_next){
                        value_next = components0[u];
                        converged = false;
                    }
                }

                // as above for the incoming edges (only directed graphs)
                if(is_directed){
                    uint64_t degree_in = view->get_degree_in(v);
                    if(degree_in > neighbours_sz){
                        COUT_DEBUG_WCC("realloc the neighbours (in): " << neighbours_sz << " -> " << degree_in);
                        neighbours_sz = degree_in;
                        neighbours = (lite_edge_t*) realloc(neighbours, sizeof(neighbours[0]) * neighbours_sz);
                        if(neighbours == nullptr) throw std::bad_alloc{};
                    }
                    view->get_nebrs_in(v, neighbours);

                    for(uint64_t i = 0; i < degree_in; i++){
                        uint64_t u = get_sid(neighbours[i]);
                        if(components0[u] < value_next){
                            value_next = components0[u];
                            converged = false;
                        }
                    }
                }

                // set the value for the next iteration
                components1[v] = value_next;
            }

            free(neighbours); neighbours = nullptr; neighbours_sz = 0;
        }

        std::swap(components0, components1);
    } while (!converged && !timeout.is_timeout());

    delete_static_view(view);


    // the two pointers may have been swapped and lost the correspondence components0 == *(ptr_components0)
    if(components0 == ptr_components0.get()){
        return ptr_components0;
    } else {
        return ptr_components1;
    }
}

void GraphOne::wcc(const char* dump2file) {
    // Init
    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();
    uint64_t N = num_vertices();

    // Run the WCC algorithm
    unique_ptr<uint32_t[]> ptr_components = graphone_execute_wcc(N, m_is_directed, timeout);
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);  }

    // Translate the vertex IDs and dump to file if required
    uint32_t* components = ptr_components.get();
    if(m_translate_vertex_ids) { // translate from GraphOne vertex ids to external vertex ids
        cuckoohash_map</* external id */ string, /* component */ uint32_t> external_ids;
        #pragma omp parallel for
        for(sid_t internal_id = 0; internal_id < N; internal_id++){
            // 1.what is the real node ID, in the external domain (e.g. user id)
            string vertex_name = g->get_typekv()->get_vertex_name(internal_id);
            if(vertex_name.empty()) continue; // this means that a mapping does not exist. It should never occur as atm we don't support vertex deletions

            // 2. retrieve the distance / weight
            uint32_t component = components[internal_id];

            // 3. make the association vertex name - score
            external_ids.insert(vertex_name, component);
        }

        if(timeout.is_timeout()){
            RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);
        }

        // store the results in the given file
        if(dump2file != nullptr){
            COUT_DEBUG("save the results to: " << dump2file);
            fstream handle(dump2file, ios_base::out);
            if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

            auto hashtable = external_ids.lock_table();

            for(const auto& keyvaluepair : hashtable){
                handle << keyvaluepair.first << " " << keyvaluepair.second << "\n";
            }

            handle.close();
        }

    } else { // without the vertex dictionary

        // store the results in the given file
        if(dump2file != nullptr){
            COUT_DEBUG("save the results to: " << dump2file);
            fstream handle(dump2file, ios_base::out);
            if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

            for(sid_t internal_id = 0; internal_id < N; internal_id++){
                handle << internal_id << " " << components[internal_id] << "\n";
            }

            handle.close();
        }
    }
}

/*****************************************************************************
 *                                                                           *
 *  CDLP                                                                      *
 *                                                                           *
 *****************************************************************************/
// same impl~ as the one done for llama
static unique_ptr<uint64_t[]> graphone_execute_cdlp(bool is_directed, bool translate_vertex_names, uint64_t num_vertices, uint64_t max_iterations, utility::TimeoutService& timer){
    auto* view = create_static_view(get_graphone_graph(), SIMPLE_MASK); // global
    unique_ptr<uint64_t[]> ptr_labels0 { new uint64_t[num_vertices] };
    unique_ptr<uint64_t[]> ptr_labels1 { new uint64_t[num_vertices] };
    uint64_t* labels0 = ptr_labels0.get(); // current labels
    uint64_t* labels1 = ptr_labels1.get(); // labels for the next iteration

    // initialisation
    if(translate_vertex_names){
        auto vertex_dictionary = g->get_typekv();

        #pragma omp parallel for shared(vertex_dictionary)
        for(uint64_t v = 0; v < num_vertices; v++){
            string str_vertex_name = vertex_dictionary->get_vertex_name(v);
            assert(!str_vertex_name.empty() && "Empty implies that the vertex does not exist, but this should not be the case because vertex "
                    "deletions are not supported atm & there should be no gaps between [0, num_vertices)");
            uint64_t vertex_name = std::stoull(str_vertex_name);
            labels0[v] = vertex_name;
        }
    } else { // dense vertex IDs, in [0, num_vertices)
        #pragma omp parallel for
        for(uint64_t v = 0; v < num_vertices; v++){
            labels0[v] = v;
        }
    }

    // algorithm pass
    bool change = true;
    uint64_t current_iteration = 0;
    while(current_iteration < max_iterations && change && !timer.is_timeout()){
        change = false; // reset the flag

        #pragma omp parallel shared(change)
        {
            lite_edge_t* neighbours = nullptr;
            uint64_t neighbours_sz = 0;

            #pragma omp for
            for(uint64_t v = 0; v < num_vertices; v++){
                unordered_map<uint64_t, uint64_t> histogram;

                // compute the histogram from both the outgoing & incoming edges. The aim is to find the number of each label
                // is shared among the neighbours of node_id
                uint64_t degree_out = view->get_degree_out(v);
                if(degree_out > neighbours_sz){
                    neighbours_sz = degree_out;
                    neighbours = (lite_edge_t*) realloc(neighbours, sizeof(neighbours[0]) * neighbours_sz);
                    if(neighbours == nullptr) throw std::bad_alloc{};
                }
                view->get_nebrs_out(v, neighbours);

                for(uint64_t i = 0; i < degree_out; i++){
                    uint64_t u = get_sid(neighbours[i]);
                    histogram[labels0[u]]++;
                }

                // cfr. Spec v0.9 pp 14 "If the graph is directed and a neighbor is reachable via both an incoming and
                // outgoing edge, its label will be counted twice"
                if(is_directed){
                    uint64_t degree_in = view->get_degree_in(v);
                    if(degree_in > neighbours_sz){
                        neighbours_sz = degree_in;
                        neighbours = (lite_edge_t*) realloc(neighbours, sizeof(neighbours[0]) * neighbours_sz);
                        if(neighbours == nullptr) throw std::bad_alloc{};
                    }
                    view->get_nebrs_in(v, neighbours);

                    for(uint64_t i = 0; i < degree_in; i++){
                        uint64_t u = get_sid(neighbours[i]);
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


            free(neighbours); neighbours = nullptr; neighbours_sz = 0;
        }

        std::swap(labels0, labels1); // next iteration
        current_iteration++;
    }

    delete_static_view(view);

    if(labels0 == ptr_labels0.get()){
        return ptr_labels0;
    } else {
        return ptr_labels1;
    }
}

void GraphOne::cdlp(uint64_t max_iterations, const char* dump2file) {
    // Init
    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();
    uint64_t N = num_vertices();

    // Run the CDLP algorithm
    unique_ptr<uint64_t[]> labels = graphone_execute_cdlp(m_is_directed, m_translate_vertex_ids, N, max_iterations, timeout);
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);  }

    // Translate the vertex IDs and dump to file if required
    if(m_translate_vertex_ids) { // translate from GraphOne vertex ids to external vertex ids
        cuckoohash_map</* external id */ string, /* component */ uint64_t> external_ids;
        #pragma omp parallel for
        for(sid_t internal_id = 0; internal_id < N; internal_id++){
            // 1.what is the real node ID, in the external domain (e.g. user id)
            string vertex_name = g->get_typekv()->get_vertex_name(internal_id);
            if(vertex_name.empty()) continue; // this means that a mapping does not exist. It should never occur as atm we don't support vertex deletions

            // 2. retrieve the distance / weight
            uint64_t component = labels[internal_id];

            // 3. make the association vertex name - score
            external_ids.insert(vertex_name, component);
        }

        if(timeout.is_timeout()){
            RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);
        }

        // store the results in the given file
        if(dump2file != nullptr){
            COUT_DEBUG("save the results to: " << dump2file);
            fstream handle(dump2file, ios_base::out);
            if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

            auto hashtable = external_ids.lock_table();

            for(const auto& keyvaluepair : hashtable){
                handle << keyvaluepair.first << " " << keyvaluepair.second << "\n";
            }

            handle.close();
        }

    } else { // without the vertex dictionary

        // store the results in the given file
        if(dump2file != nullptr){
            COUT_DEBUG("save the results to: " << dump2file);
            fstream handle(dump2file, ios_base::out);
            if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

            for(sid_t internal_id = 0; internal_id < N; internal_id++){
                handle << internal_id << " " << labels[internal_id] << "\n";
            }

            handle.close();
        }
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

// loosely based on the impl~ made for Stinger
static unique_ptr<double[]> graphone_execute_lcc(uint64_t num_vertices, bool is_directed, utility::TimeoutService& timer){
    auto* view = create_static_view(get_graphone_graph(), SIMPLE_MASK); // global
    unique_ptr<double[]> ptr_lcc { new double[num_vertices] };
    double* lcc = ptr_lcc.get();

    #pragma omp parallel
    {

        // to iterate over the sequence of edges
        lite_edge_t* edges1 = nullptr;
        uint64_t edges1_sz = 0;
        lite_edge_t* edges2 = nullptr;
        uint64_t edges2_sz = 0;

        #pragma omp for
        for(uint64_t v = 0; v < num_vertices; v++){
            COUT_DEBUG_LCC("> Node " << v);
            if(timer.is_timeout()) continue; // exhausted the budget of available time
            lcc[v] = 0.0;
            uint64_t num_triangles = 0; // number of triangles found so far for the node v

            // Cfr. Spec v.0.9.0 pp. 15: "If the number of neighbors of a vertex is less than two, its coefficient is defined as zero"
            uint64_t v_degree_out = view->get_degree_out(v);
            uint64_t v_degree_in = is_directed ? view->get_degree_in(v) : 0;
            uint64_t v_degree_ub = v_degree_in + v_degree_out; // upper bound for directed graphs, exact degree for those undirected
            if(v_degree_ub < 2) continue;

            // Build the list of neighbours of v
            unordered_set<uint64_t> neighbours;
            uint64_t v_degree = v_degree_out; // actual degree, considering the intersection of N_in(v) and N_out(v)
            if(edges1_sz < v_degree_ub){
                COUT_DEBUG_LCC("* realloc (1): " << edges1_sz << " -> " << v_degree_ub);
                edges1_sz = v_degree_ub;
                edges1 = reinterpret_cast<decltype(edges1)>( realloc(edges1, sizeof(edges1[0]) * edges1_sz) );
                if(edges1 == nullptr) throw std::bad_alloc{};
            }


            // Outgoing edges
            view->get_nebrs_out(v, edges1);
            for(uint64_t i = 0; i < v_degree_out; i++){
                uint64_t u = get_sid(edges1[i]);
                neighbours.insert(u);
            }

            // Incoming edges (only directed graphs)
            if(is_directed){
                if(edges2_sz < v_degree_in){
                    COUT_DEBUG_LCC("* realloc (2.1): " << edges2_sz << " -> " << v_degree_in);
                    edges2_sz = v_degree_in;
                    edges2 = reinterpret_cast<decltype(edges2)>( realloc(edges2, sizeof(edges2[0]) * edges2_sz) );
                    if(edges2 == nullptr) throw std::bad_alloc{};
                }

                view->get_nebrs_in(v, edges2);
                for(uint64_t i = 0; i < v_degree_in; i++){
                    uint64_t u = get_sid(edges2[i]);
                    auto result = neighbours.insert(u);
                    if(result.second){ // the element was actually inserted
                        edges1[v_degree++] = edges2[i];
                    }
                }
            }

            // Now we know is the actual degree of v, perform the proper check for directed graphs
            if(is_directed && v_degree < 2) continue;

            // again, visit all neighbours of v
            // for directed graphs, edges1 contains the intersection of both the incoming and the outgoing edges
            for(uint64_t i = 0; i < v_degree; i++){
                uint64_t u = get_sid(edges1[i]);
                COUT_DEBUG_LCC("[" << i << "/" << v_degree << "] neighbour: " << u);
                assert(neighbours.count(u) == 1 && "The set `neighbours' should contain all neighbours of v");

                uint64_t u_degree_out = view->get_degree_out(u);
                if(edges2_sz < u_degree_out){
                    COUT_DEBUG_LCC("* realloc (2.2): " << edges2_sz << " -> " << u_degree_out);
                    edges2_sz = u_degree_out;
                    edges2 = reinterpret_cast<decltype(edges2)>( realloc(edges2, sizeof(edges2[0]) * edges2_sz) );
                    if(edges2 == nullptr) throw std::bad_alloc{};
                }

                // For the Graphalytics spec v 0.9.0, only consider the outgoing edges for the neighbours u
                view->get_nebrs_out(u, edges2);
                for(uint64_t j = 0; j < u_degree_out; j++){
                    uint64_t w = get_sid(edges2[j]);
                    COUT_DEBUG_LCC("---> [" << j << "/" << u_degree_out << "] neighbour: " << w);
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

        free(edges1); edges1 = nullptr; edges1_sz = 0;
        free(edges2); edges2 = nullptr; edges2_sz = 0;
    }

    delete_static_view(view);

    return ptr_lcc;
}

void GraphOne::lcc(const char* dump2file) {
    // Init
    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();
    uint64_t N = num_vertices();

    // Run the LCC algorithm
    unique_ptr<double[]> scores = graphone_execute_lcc(N, is_directed(), timeout);
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);  }

    // Translate the vertex IDs and dump to file if required
    if(m_translate_vertex_ids) { // translate from GraphOne vertex ids to external vertex ids
        cuckoohash_map</* external id */ string, /* score */ double> external_ids;
        #pragma omp parallel for
        for(sid_t internal_id = 0; internal_id < N; internal_id++){
            // 1.what is the real node ID, in the external domain (e.g. user id)
            string vertex_name = g->get_typekv()->get_vertex_name(internal_id);
            if(vertex_name.empty()) continue; // this means that a mapping does not exist. It should never occur as atm we don't support vertex deletions

            // 2. retrieve the distance / weight
            double score = scores[internal_id];

            // 3. make the association vertex name - score
            external_ids.insert(vertex_name, score);
        }

        if(timeout.is_timeout()){
            RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);
        }

        // store the results in the given file
        if(dump2file != nullptr){
            COUT_DEBUG("save the results to: " << dump2file);
            fstream handle(dump2file, ios_base::out);
            if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

            auto hashtable = external_ids.lock_table();

            for(const auto& keyvaluepair : hashtable){
                handle << keyvaluepair.first << " " << keyvaluepair.second << "\n";
            }

            handle.close();
        }

    } else { // without the vertex dictionary

        // store the results in the given file
        if(dump2file != nullptr){
            COUT_DEBUG("save the results to: " << dump2file);
            fstream handle(dump2file, ios_base::out);
            if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

            for(sid_t internal_id = 0; internal_id < N; internal_id++){
                handle << internal_id << " " << scores[internal_id] << "\n";
            }

            handle.close();
        }
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

namespace { // anonymous

/*
GAP Benchmark Suite
Class:  pvector
Author: Scott Beamer

Vector class with ability to not initialize or do initialize in parallel
 - std::vector (when resizing) will always initialize, and does it serially
 - When pvector is resized, new elements are uninitialized
 - Resizing is not thread-safe
*/

template <typename T_>
class pvector {
 public:
  typedef T_* iterator;

  pvector() : start_(nullptr), end_size_(nullptr), end_capacity_(nullptr) {}

  explicit pvector(size_t num_elements) {
    start_ = new T_[num_elements];
    end_size_ = start_ + num_elements;
    end_capacity_ = end_size_;
  }

  pvector(size_t num_elements, T_ init_val) : pvector(num_elements) {
    fill(init_val);
  }

  pvector(iterator copy_begin, iterator copy_end)
      : pvector(copy_end - copy_begin) {
    #pragma omp parallel for
    for (size_t i=0; i < capacity(); i++)
      start_[i] = copy_begin[i];
  }

  // don't want this to be copied, too much data to move
  pvector(const pvector &other) = delete;

  // prefer move because too much data to copy
  pvector(pvector &&other)
      : start_(other.start_), end_size_(other.end_size_),
        end_capacity_(other.end_capacity_) {
    other.start_ = nullptr;
    other.end_size_ = nullptr;
    other.end_capacity_ = nullptr;
  }

  // want move assignment
  pvector& operator= (pvector &&other) {
    start_ = other.start_;
    end_size_ = other.end_size_;
    end_capacity_ = other.end_capacity_;
    other.start_ = nullptr;
    other.end_size_ = nullptr;
    other.end_capacity_ = nullptr;
    return *this;
  }

  ~pvector() {
    if (start_ != nullptr)
      delete[] start_;
  }

  // not thread-safe
  void reserve(size_t num_elements) {
    if (num_elements > capacity()) {
      T_ *new_range = new T_[num_elements];
      #pragma omp parallel for
      for (size_t i=0; i < size(); i++)
        new_range[i] = start_[i];
      end_size_ = new_range + size();
      delete[] start_;
      start_ = new_range;
      end_capacity_ = start_ + num_elements;
    }
  }

  bool empty() {
    return end_size_ == start_;
  }

  void clear() {
    end_size_ = start_;
  }

  void resize(size_t num_elements) {
    reserve(num_elements);
    end_size_ = start_ + num_elements;
  }

  T_& operator[](size_t n) {
    return start_[n];
  }

  const T_& operator[](size_t n) const {
    return start_[n];
  }

  void push_back(T_ val) {
    if (size() == capacity()) {
      size_t new_size = capacity() == 0 ? 1 : capacity() * growth_factor;
      reserve(new_size);
    }
    *end_size_ = val;
    end_size_++;
  }

  void fill(T_ init_val) {
    #pragma omp parallel for
    for (T_* ptr=start_; ptr < end_size_; ptr++)
      *ptr = init_val;
  }

  size_t capacity() const {
    return end_capacity_ - start_;
  }

  size_t size() const {
    return end_size_ - start_;
  }

  iterator begin() const {
    return start_;
  }

  iterator end() const {
    return end_size_;
  }

  T_* data() const {
    return start_;
  }

  void swap(pvector &other) {
    std::swap(start_, other.start_);
    std::swap(end_size_, other.end_size_);
    std::swap(end_capacity_, other.end_capacity_);
  }


 private:
  T_* start_;
  T_* end_size_;
  T_* end_capacity_;
  static const size_t growth_factor = 2;
};

/*
GAP Benchmark Suite
File:   Platform Atomics
Author: Scott Beamer

Wrappers for compiler intrinsics for atomic memory operations (AMOs)
 - If not using OpenMP (serial), provides serial fallbacks
*/

#if defined _OPENMP

  #if defined __GNUC__

    // gcc/clang/icc instrinsics

    template<typename T, typename U>
    T fetch_and_add(T &x, U inc) {
      return __sync_fetch_and_add(&x, inc);
    }

    template<typename T>
    bool compare_and_swap(T &x, const T &old_val, const T &new_val) {
      return __sync_bool_compare_and_swap(&x, old_val, new_val);
    }

    template<>
    [[maybe_unused]] bool compare_and_swap(float &x, const float &old_val, const float &new_val) {
      return __sync_bool_compare_and_swap(reinterpret_cast<uint32_t*>(&x),
                                          reinterpret_cast<const uint32_t&>(old_val),
                                          reinterpret_cast<const uint32_t&>(new_val));
    }

    template<>
    [[maybe_unused]] bool compare_and_swap(double &x, const double &old_val, const double &new_val) {
      return __sync_bool_compare_and_swap(reinterpret_cast<uint64_t*>(&x),
                                          reinterpret_cast<const uint64_t&>(old_val),
                                          reinterpret_cast<const uint64_t&>(new_val));
    }

  #elif __SUNPRO_CC

    // sunCC (solaris sun studio) intrinsics
    // less general, only work for int32_t, int64_t, uint32_t, uint64_t
    // http://docs.oracle.com/cd/E19253-01/816-5168/6mbb3hr06/index.html

    #include <atomic.h>
    #include <cinttypes>

    int32_t fetch_and_add(int32_t &x, int32_t inc) {
      return atomic_add_32_nv((volatile uint32_t*) &x, inc) - inc;
    }

    int64_t fetch_and_add(int64_t &x, int64_t inc) {
      return atomic_add_64_nv((volatile uint64_t*) &x, inc) - inc;
    }

    uint32_t fetch_and_add(uint32_t &x, uint32_t inc) {
      return atomic_add_32_nv((volatile uint32_t*) &x, inc) - inc;
    }

    uint64_t fetch_and_add(uint64_t &x, uint64_t inc) {
      return atomic_add_64_nv((volatile uint64_t*) &x, inc) - inc;
    }

    bool compare_and_swap(int32_t &x, const int32_t &old_val, const int32_t &new_val) {
      return old_val == atomic_cas_32((volatile uint32_t*) &x, old_val, new_val);
    }

    bool compare_and_swap(int64_t &x, const int64_t &old_val, const int64_t &new_val) {
      return old_val == atomic_cas_64((volatile uint64_t*) &x, old_val, new_val);
    }

    bool compare_and_swap(uint32_t &x, const uint32_t &old_val, const uint32_t &new_val) {
      return old_val == atomic_cas_32((volatile uint32_t*) &x, old_val, new_val);
    }

    bool compare_and_swap(uint64_t &x, const uint64_t &old_val, const uint64_t &new_val) {
      return old_val == atomic_cas_64((volatile uint64_t*) &x, old_val, new_val);
    }

    bool compare_and_swap(float &x, const float &old_val, const float &new_val) {
      return old_val == atomic_cas_32((volatile uint32_t*) &x,
                                      (const volatile uint32_t&) old_val,
                                      (const volatile uint32_t&) new_val);
    }

    bool compare_and_swap(double &x, const double &old_val, const double &new_val) {
      return old_val == atomic_cas_64((volatile uint64_t*) &x,
                                      (const volatile uint64_t&) old_val,
                                      (const volatile uint64_t&) new_val);
    }

  #else   // defined __GNUC__ __SUNPRO_CC

    #error No atomics available for this compiler but using OpenMP

  #endif  // else defined __GNUC__ __SUNPRO_CC

#else   // defined _OPENMP

  // serial fallbacks

  template<typename T, typename U>
  T fetch_and_add(T &x, U inc) {
    T orig_val = x;
    x += inc;
    return orig_val;
  }

  template<typename T>
  bool compare_and_swap(T &x, const T &old_val, const T &new_val) {
    if (x == old_val) {
      x = new_val;
      return true;
    }
    return false;
  }

#endif  // else defined _OPENMP

} // anonymous namespace

using NodeID = uint64_t;
using WeightT = double;
static const size_t kMaxBin = numeric_limits<size_t>::max()/2;

static pvector<WeightT> graphone_execute_sssp(uint64_t num_vertices, uint64_t source, double delta, utility::TimeoutService& timer){
    auto* view = create_static_view(get_graphone_graph(), SIMPLE_MASK); // global

    // Total number of edges in the view
    uint64_t num_edges = 0;
    #pragma omp parallel for reduction(+:num_edges)
    for(uint64_t v = 0; v < num_vertices; v++){
        num_edges += view->get_degree_out(v);
    }

    // Init
    pvector<WeightT> dist(num_vertices, numeric_limits<WeightT>::infinity());
    dist[source] = 0;
    pvector<NodeID> frontier(num_edges);
    // two element arrays for double buffering curr=iter&1, next=(iter+1)&1
    size_t shared_indexes[2] = {0, kMaxBin};
    size_t frontier_tails[2] = {1, 0};
    frontier[0] = source;

    #pragma omp parallel
    {
        vector<vector<NodeID> > local_bins(0);
        size_t iter = 0;
        lite_edge_t* neighbours = nullptr;
        uint64_t neighbours_sz = 0;

        while (shared_indexes[iter&1] != kMaxBin) {
            size_t &curr_bin_index = shared_indexes[iter&1];
            size_t &next_bin_index = shared_indexes[(iter+1)&1];
            size_t &curr_frontier_tail = frontier_tails[iter&1];
            size_t &next_frontier_tail = frontier_tails[(iter+1)&1];
            #pragma omp for nowait schedule(dynamic, 64)
            for (size_t i=0; i < curr_frontier_tail; i++) {
                NodeID u = frontier[i];
                if (dist[u] >= delta * static_cast<WeightT>(curr_bin_index)) {

                    uint64_t u_degree = view->get_degree_out(u);
                    if(u_degree > neighbours_sz){
                        neighbours_sz = u_degree;
                        neighbours = reinterpret_cast<decltype(neighbours)>( realloc(neighbours, sizeof(neighbours[0]) * neighbours_sz));
                        if(neighbours == nullptr) throw std::bad_alloc{};
                    }

                    view->get_nebrs_out(u, neighbours);
                    for(uint64_t i = 0; i < u_degree; i++){
                        uint64_t v = get_sid(neighbours[i]);
                        double w = neighbours[i].second.value_double;

                        WeightT old_dist = dist[v];
                        WeightT new_dist = dist[u] + w;
                        if (new_dist < old_dist) {
                            bool changed_dist = true;
                            while (!compare_and_swap(dist[v], old_dist, new_dist)) {
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
                size_t copy_start = fetch_and_add(next_frontier_tail, local_bins[next_bin_index].size());
                copy(local_bins[next_bin_index].begin(), local_bins[next_bin_index].end(), frontier.data() + copy_start);
                local_bins[next_bin_index].resize(0);
            }

            iter++;
        #pragma omp barrier
        }

        free(neighbours); neighbours = nullptr; neighbours_sz = 0;

#if defined(DEBUG)
        #pragma omp single
        COUT_DEBUG("took " << iter << " iterations");
#endif
    }

    delete_static_view(view);
    return dist;
}

void GraphOne::sssp(uint64_t source_vertex_id, const char* dump2file) {
    // Init
    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();
    uint64_t root = vtx_ext2int(source_vertex_id); // from a graphalytics vertex ID to the dense vertex ID used internally
    uint64_t N = num_vertices();

    // Run the SSSP algorithm
    double delta = 2.0; // same value used in the GAPBS, at least for most graphs
    auto distances = graphone_execute_sssp(N, root, delta, timeout);
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);  }

    // Translate the vertex IDs and dump to file if required
    if(m_translate_vertex_ids) { // translate from GraphOne vertex ids to external vertex ids
        cuckoohash_map</* external id */ string, /* distance */ double> external_ids;
#pragma omp parallel for
        for(sid_t internal_id = 0; internal_id < N; internal_id++){
            // 1.what is the real node ID, in the external domain (e.g. user id)
            string vertex_name = g->get_typekv()->get_vertex_name(internal_id);
            if(vertex_name.empty()) continue; // this means that a mapping does not exist. It should never occur as atm we don't support vertex deletions

            // 2. retrieve the distance / weight
            double distance = distances[internal_id];

            // 3. make the association vertex name - distance
            external_ids.insert(vertex_name, distance);
        }

        if(timeout.is_timeout()){
            RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);
        }

        // store the results in the given file
        if(dump2file != nullptr){
            COUT_DEBUG("save the results to: " << dump2file);
            fstream handle(dump2file, ios_base::out);
            if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

            auto hashtable = external_ids.lock_table();

            for(const auto& keyvaluepair : hashtable){
                handle << keyvaluepair.first << " " << keyvaluepair.second << "\n";
            }

            handle.close();
        }

    } else { // without the vertex dictionary

        // store the results in the given file
        if(dump2file != nullptr){
            COUT_DEBUG("save the results to: " << dump2file);
            fstream handle(dump2file, ios_base::out);
            if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");


            for(sid_t internal_id = 0; internal_id < N; internal_id++){
                handle << internal_id << " " << distances[internal_id] << "\n";
            }

            handle.close();
        }
    }
}

} // namespace


