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

#include "llama_class.hpp"
#include "llama_internal.hpp"

#include <cmath>
#include <iostream>
#include <mutex>
#include <shared_mutex> // shared_lock

#include "common/time.hpp"

using namespace common;
using namespace std;

// External counters used by LLAMA
#if defined(LL_COUNTERS)
std::atomic<size_t> g_iter_begin;
std::atomic<size_t> g_iter_descend;
std::atomic<size_t> g_iter_next;
#endif

// External counters, to profile ll_writable_graph#add_edge_if_not_exists and #build
// defined in llama_internal.hpp
#if defined(LL_PROFILE_UPDATES)
common::Timer<true> g_llama_total_time;
std::atomic<uint64_t> g_llama_add_edge_check_nanosecs = 0;
std::atomic<uint64_t> g_llama_add_edge_total_nanosecs = 0;
std::atomic<uint64_t> g_llama_build_nanosecs = 0;
void llama_add_edge_print_stats();
#endif

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG
namespace gfe { extern mutex _log_mutex [[maybe_unused]]; }
#define COUT_DEBUG_FORCE(msg) { std::scoped_lock<std::mutex> lock{gfe::_log_mutex}; std::cout << "[LLAMA::" << __FUNCTION__ << "] " << msg << std::endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif

/*****************************************************************************
 *                                                                           *
 *  Helpers                                                                  *
 *                                                                           *
 *****************************************************************************/
// Get the weight associated to an _outgoing_ edge_id
template<typename Graph>
static double get_out_edge_weight(Graph& graph, edge_t edge_id){
    return reinterpret_cast<ll_mlcsr_edge_property<double>*>(graph.get_edge_property_64(g_llama_property_weights))->get(edge_id);
}

// Retrieve the weight associated to the given edge_id, assuming that it is an incoming edge
// In the read-only store, the weight is only a property of LLAMA's outgoing edges. Retrieve the corresponding edge b -> a from a -> b
static double get_in_edge_weight(ll_mlcsr_ro_graph& graph, edge_t edge_id){
    return get_out_edge_weight(graph,  graph.in_to_out(edge_id) );
}

// Retrieve the weight associated to the given edge_id, assuming that it is an incoming edge
//  In the delta store, edge_ids are the same for both incoming & outgoing edges. We can discriminate between read-only and delta
// edges because delta store edges have a negative edge_id
static double get_in_edge_weight(ll_writable_graph& graph, edge_t edge_id){
    if(edge_id < 0){ // write store
        return get_out_edge_weight(graph, edge_id);
    } else { // read store
        return get_in_edge_weight(graph.ro_graph(), edge_id);
    }
}

#if defined(LL_PROFILE_UPDATES)
void llama_add_edge_print_stats() {
    uint64_t total_time_nanosecs = g_llama_total_time.nanoseconds();

    double oh_build = g_llama_build_nanosecs * 100.0 / total_time_nanosecs;

    COUT_DEBUG_FORCE("Total time in updates: " << common::time::to_string(chrono::nanoseconds(total_time_nanosecs)) );
    COUT_DEBUG_FORCE("Cumulative time spent in #build: " << common::time::to_string(chrono::nanoseconds(g_llama_build_nanosecs) ) << ", percentage: " << oh_build << "%");

    double oh_add_edge = g_llama_add_edge_total_nanosecs * 100.0 / total_time_nanosecs;
    COUT_DEBUG_FORCE("Cumulative time spent in #add_edge_if_not_exists: " << common::time::to_string(chrono::nanoseconds( g_llama_add_edge_total_nanosecs ) ) << ", percentage: " << oh_add_edge << "%");

    double oh_add_edge_check = g_llama_add_edge_check_nanosecs * 100.0 / g_llama_add_edge_total_nanosecs;
    COUT_DEBUG_FORCE("Cumulative time spent checking edge existence in #add_edge_if_not_exists: " << common::time::to_string(chrono::nanoseconds( g_llama_add_edge_check_nanosecs ) ) << ", overhead: " << oh_add_edge_check << "%");
}
#endif

/*****************************************************************************
 *                                                                           *
 *  Init                                                                     *
 *                                                                           *
 *****************************************************************************/

namespace gfe::library {

LLAMAClass::LLAMAClass(bool is_directed) : m_is_directed(is_directed) {
    m_db = new ll_database();

    auto& csr = m_db->graph()->ro_graph();
    csr.create_uninitialized_node_property_64(g_llama_property_names, LL_T_INT64); // translate a logical id (e.g. node_id = 4) into an external node id (e.g. user_id = 21084718954)
    csr.create_uninitialized_edge_property_64(g_llama_property_weights, LL_T_DOUBLE);

#if !defined(LLAMA_HASHMAP_WITH_TBB)
    m_vmap_locks = new SpinLock[m_num_vmap_locks];
#endif
}


LLAMAClass::~LLAMAClass(){
#if defined(LL_PROFILE_UPDATES)
    llama_add_edge_print_stats();
#endif

#if !defined(LLAMA_HASHMAP_WITH_TBB)
    delete[] m_vmap_locks; m_vmap_locks = nullptr;
#endif
    delete m_db; m_db = nullptr;
}

/*****************************************************************************
 *                                                                           *
 *  Properties                                                               *
 *                                                                           *
 *****************************************************************************/
uint64_t LLAMAClass::num_edges() const {
    return m_num_edges;
}

uint64_t LLAMAClass::num_vertices() const {
#if defined(LLAMA_HASHMAP_WITH_TBB)
    return m_vmap.size();
#else
    return m_vmap_read_only.size();
#endif
}

uint64_t LLAMAClass::num_levels() const{
    shared_lock<shared_mutex_t> cplock(m_lock_checkpoint); // forbid any checkpoint now
    return m_db->graph()->ro_graph().num_levels();
}

bool LLAMAClass::is_directed() const {
    return m_is_directed;
}

bool LLAMAClass::has_vertex(uint64_t vertex_id) const {
#if defined(LLAMA_HASHMAP_WITH_TBB)
    decltype(m_vmap)::const_accessor accessor;
    return m_vmap.find(accessor, vertex_id);
#else
    scoped_lock<SpinLock> vlock(m_vmap_locks[vertex_id % m_num_vmap_locks]);

    try {
        vmap_write_store_find(vertex_id); // exception -> vertex not present
        return true;
    } catch (std::out_of_range& e){
        return false;
    }
#endif
}


double LLAMAClass::get_weight(uint64_t ext_source_id, uint64_t ext_destination_id) const {
    auto nan = numeric_limits<double>::quiet_NaN();


#if defined(LLAMA_HASHMAP_WITH_TBB)
    decltype(m_vmap)::const_accessor accessor1, accessor2; // shared lock on the dictionary
    if(!m_vmap.find(accessor1, ext_source_id)){ return nan; }
    if(!m_vmap.find(accessor2, ext_destination_id)) { return nan; }

    int64_t source = accessor1->second;
    int64_t destination = accessor2->second;

    if(!m_is_directed && source > destination){ // for undirected graphs, ensure src_id < dst_id
        std::swap(source, destination);
    }

    edge_t edge_id = m_db->graph()->find(source, destination);
    if(edge_id == LL_NIL_EDGE) return nan; // the edge does not exist

    return get_out_edge_weight(* (m_db->graph()), edge_id);

#else
    uint64_t int_source_id, int_destination_id;
    try {
        uint64_t vtx_lock_1 = min(ext_source_id % m_num_vmap_locks, ext_destination_id % m_num_vmap_locks);
        uint64_t vtx_lock_2 = max(ext_source_id % m_num_vmap_locks, ext_destination_id % m_num_vmap_locks);
        scoped_lock<SpinLock> vlock1(m_vmap_locks[vtx_lock_1]);
        scoped_lock<SpinLock> vlock2(m_vmap_locks[vtx_lock_2]);

        int_source_id = vmap_write_store_find(ext_source_id); // throws std::out_of_range if ext_source_id is not present
        int_destination_id = vmap_write_store_find(ext_destination_id); // as above

        if(!m_is_directed){ // for undirected graphs, ensure src_id < dst_id
            if(int_source_id > int_destination_id){
                std::swap(int_source_id, int_destination_id);
            }
        }

        edge_t edge_id = m_db->graph()->find(int_source_id, int_destination_id);
        if(edge_id == LL_NIL_EDGE) return nan; // the edge does not exist

        return get_out_edge_weight(* (m_db->graph()), edge_id);
    } catch ( std::out_of_range& e){
        return nan; // either ext_source_id or ext_destination_id
    }
#endif
}

void LLAMAClass::set_timeout(uint64_t seconds) {
    m_timeout = chrono::seconds{ seconds };
}


#if !defined(LLAMA_HASHMAP_WITH_TBB)
int64_t LLAMAClass::vmap_write_store_find(uint64_t external_vertex_id) const {
    int64_t logical_id = -1;

    if( m_vmap_updated.find(external_vertex_id, logical_id) ){ // found
        return logical_id;
    }

    if( m_vmap_removed.contains(external_vertex_id) ){ // vertex explicitly removed
        throw std::out_of_range("vertex deleted");
    }

    if( m_vmap_read_only.find(external_vertex_id, logical_id) ){ // found
        return logical_id;
    }

    throw std::out_of_range("the mapping is not present");
}

bool LLAMAClass::vmap_write_store_contains(uint64_t external_vertex_id) const {
    try {
        vmap_write_store_find(external_vertex_id); // throws std::out_of_range if not present
        return true;
    } catch (std::out_of_range& e){
        return false;
    }
}
#endif

int64_t LLAMAClass::get_internal_vertex_id(uint64_t external_vertex_id) const {
#if defined(LLAMA_HASHMAP_WITH_TBB)
    decltype(m_vmap)::const_accessor accessor;
    if ( m_vmap.find(accessor, external_vertex_id ) ){
        return accessor->second;
    } else {
        ERROR("The given vertex does not exist: " << external_vertex_id);
    }
#else
    int64_t llama_source_vertex_id;
    if (! m_vmap_read_only.find(external_vertex_id, llama_source_vertex_id) ){ // side effect, it assigns llama_source_vertex_id
        ERROR("The given vertex does not exist: " << external_vertex_id);
    }
#endif
}

uint64_t LLAMAClass::get_read_store_outdegree(ll_mlcsr_ro_graph& snapshot, int64_t llama_vertex_id) const{
    if(m_is_directed){
        return snapshot.out_degree(llama_vertex_id);
    } else {
        return snapshot.out_degree(llama_vertex_id) + snapshot.in_degree(llama_vertex_id);
    }
}

uint64_t LLAMAClass::get_write_store_outdegree(int64_t llama_vertex_id) const{
    ll_writable_graph* graph = m_db->graph();

    if(m_is_directed){
        return graph->out_degree(llama_vertex_id);
    } else {
        return graph->out_degree(llama_vertex_id) + graph->in_degree(llama_vertex_id);
    }
}

ll_mlcsr_ro_graph LLAMAClass::get_snapshot() const {
    if(num_levels() == 0)
        ERROR("There are no levels/deltas/snapshots available. Create them with #build()");

    return ll_mlcsr_ro_graph{ &(m_db->graph()->ro_graph()) , static_cast<int>(num_levels()) -1 };
}


/*****************************************************************************
 *                                                                           *
 *  Updates                                                                  *
 *                                                                           *
 *****************************************************************************/

bool LLAMAClass::add_vertex(uint64_t vertex_id_ext){
    shared_lock<shared_mutex_t> cplock(m_lock_checkpoint); // forbid any checkpoint now
    COUT_DEBUG("vertex_id: " << vertex_id_ext);

#if defined(LLAMA_HASHMAP_WITH_TBB)
    decltype(m_vmap)::accessor accessor; // xlock
    bool inserted = m_vmap.insert(accessor, vertex_id_ext);
    if ( inserted ){
        node_t node_id = m_db->graph()->add_node();
        m_db->graph()->get_node_property_64(g_llama_property_names)->set(node_id, vertex_id_ext); // register the mapping internal -> external

        accessor->second = node_id;
    }

    return inserted;
#else
    auto& mutex = m_vmap_locks[vertex_id_ext % m_num_vmap_locks];
    mutex.lock();
    if(!vmap_write_store_contains(vertex_id_ext)){
        mutex.unlock();

        node_t node_id = m_db->graph()->add_node();
        m_db->graph()->get_node_property_64(g_llama_property_names)->set(node_id, vertex_id_ext); // register the mapping internal -> external

        /**
         * Here there may be a race condition, where both T1 and T2 invoke m_db->graph()->add_node();
         * However we roll back the effects of the second thread in #m_vmap_updated.upsert. The lambda is invoked
         * if a key is already present in the hash table (the node_id set by T1).
         */
        bool inserted = true;
        m_vmap_updated.upsert(vertex_id_ext, [this, node_id, &inserted](int64_t& previous_node_id /* ignore */){
            // roll back
            m_db->graph()->delete_node(node_id);

            inserted = false;
        }, node_id);

        return inserted;
    } else {
        // the mapping already exists
        mutex.unlock();
        return false;
    }
#endif
}

int64_t LLAMAClass::get_or_create_vertex_id(uint64_t vertex_id_ext) {
#if defined(LLAMA_HASHMAP_WITH_TBB)
    decltype(m_vmap)::accessor accessor; // xlock
    bool inserted = m_vmap.insert(accessor, vertex_id_ext);
    if ( inserted ){
        node_t node_id = m_db->graph()->add_node();
        m_db->graph()->get_node_property_64(g_llama_property_names)->set(node_id, vertex_id_ext); // register the mapping internal -> external

        accessor->second = node_id;
        return node_id;
    } else { // already exists
        return accessor->second;
    }
#else
    auto& mutex = m_vmap_locks[vertex_id_ext % m_num_vmap_locks];
    mutex.lock();
    try {
        int64_t logical_id = vmap_write_store_find(vertex_id_ext);
        mutex.unlock();

        return logical_id;

    } catch( std::out_of_range& ){ // the vertex does not exist
        mutex.unlock();

        node_t node_id = m_db->graph()->add_node();
        m_db->graph()->get_node_property_64(g_llama_property_names)->set(node_id, vertex_id_ext); // register the mapping internal -> external

        /**
         * Here there may be a race condition, where both T1 and T2 invoke m_db->graph()->add_node();
         * However we roll back the effects of the second thread in #m_vmap_updated.upsert. The lambda is invoked
         * if a key is already present in the hash table (the node_id set by T1).
         */
        m_vmap_updated.upsert(vertex_id_ext, [this, &node_id](int64_t previous_node_id){
            // the key `vertex_id_ext' already exists ...

            // roll back
            m_db->graph()->delete_node(node_id);

            // retrieve the current node_id
            node_id = previous_node_id;
        }, node_id);

        return node_id;
    }
#endif
}

bool LLAMAClass::remove_vertex(uint64_t vertex_id_ext){
    shared_lock<shared_mutex_t> cplock(m_lock_checkpoint); // forbid any checkpoint now
    COUT_DEBUG("vertex_id: " << vertex_id_ext);

#if defined(LLAMA_HASHMAP_WITH_TBB)
    decltype(m_vmap)::accessor accessor; // xlock
    bool found = m_vmap.find(accessor, vertex_id_ext);
    if(found){
        int64_t llama_vertex_id = accessor->second;
        assert(llama_vertex_id != -1);
        m_db->graph()->delete_node(llama_vertex_id);

        m_vmap.erase(accessor);
    }
    return found;
#else
    int64_t llama_vertex_id = -1;
    bool is_removed = false;

    auto& mutex = m_vmap_locks[vertex_id_ext % m_num_vmap_locks];
    mutex.lock();
    if(m_vmap_read_only.find(vertex_id_ext, llama_vertex_id)){
        m_vmap_removed.insert(vertex_id_ext, true);
        is_removed = true;
    }

    is_removed |= m_vmap_updated.erase_fn(vertex_id_ext, [&llama_vertex_id](int64_t& mapping){
        llama_vertex_id = mapping;
        return true;
    });


    mutex.unlock();

    if(is_removed){
        assert(llama_vertex_id != -1);
        m_db->graph()->delete_node(llama_vertex_id);
    }

    return is_removed;
#endif
}

bool LLAMAClass::add_edge(graph::WeightedEdge e){
    shared_lock<shared_mutex_t> cplock(m_lock_checkpoint); // forbid any checkpoint now
    COUT_DEBUG("edge: " << e);

    node_t llama_source_id { -1 }, llama_destination_id { -1 };

#if defined(LLAMA_HASHMAP_WITH_TBB)
    decltype(m_vmap)::const_accessor accessor1, accessor2; // shared lock on the dictionary
    if(!m_vmap.find(accessor1, e.source())){ return false; }
    if(!m_vmap.find(accessor2, e.destination())) { return false; }

    llama_source_id = accessor1->second;
    llama_destination_id = accessor2->second;

#else

    try {
        scoped_lock<SpinLock> vlock(m_vmap_locks[e.m_source % m_num_vmap_locks]);
        llama_source_id = vmap_write_store_find(e.m_source); // throws std::out_of_range if ext_source_id is not present
    } catch( std::out_of_range& e ){
        return false; // the source does not exist
    }

    try {
        scoped_lock<SpinLock> vlock(m_vmap_locks[e.m_destination % m_num_vmap_locks]);
        llama_destination_id = vmap_write_store_find(e.m_destination); // throws std::out_of_range if ext_source_id is not present
    } catch( std::out_of_range& e ){
        return false; // the destination does not exist
    }

#endif

    return add_edge0(llama_source_id, llama_destination_id, e.weight());
}

bool LLAMAClass::add_edge_v2(graph::WeightedEdge e){
    shared_lock<shared_mutex_t> cplock(m_lock_checkpoint); // forbid any checkpoint now
    COUT_DEBUG("edge: " << e);

    node_t llama_source_id = get_or_create_vertex_id(e.source());
    node_t llama_destination_id = get_or_create_vertex_id(e.destination());

    return add_edge0(llama_source_id, llama_destination_id, e.weight());
}

bool LLAMAClass::add_edge0(int64_t llama_source_id, int64_t llama_destination_id, double weight){
    if(!m_is_directed){ // for undirected graphs, ensure llama_source_id < llama_destination_id
        if(llama_source_id > llama_destination_id){
            std::swap(llama_source_id, llama_destination_id);
        }
    }

    edge_t edge_id;
    bool inserted = m_db->graph()->add_edge_if_not_exists(llama_source_id, llama_destination_id, &edge_id);
    assert( m_db->graph()->find(llama_source_id, llama_destination_id) == edge_id );

    // thread unsafe, this should really still be under the same latch of add_edge_if_not_exists
    if(inserted){
        m_db->graph()->get_edge_property_64(g_llama_property_weights)->set(edge_id, *reinterpret_cast<uint64_t*>(&(weight)));
    }

    return inserted;
}

bool LLAMAClass::remove_edge(graph::Edge e){
    shared_lock<shared_mutex_t> cplock(m_lock_checkpoint); // forbid any checkpoint now
    COUT_DEBUG("edge: " << e);

    node_t llama_source_id { -1 }, llama_destination_id { -1 };

#if defined(LLAMA_HASHMAP_WITH_TBB)
    decltype(m_vmap)::const_accessor accessor1, accessor2; // shared lock on the dictionary
    if(!m_vmap.find(accessor1, e.source())){ return false; }
    if(!m_vmap.find(accessor2, e.destination())) { return false; }

    llama_source_id = accessor1->second;
    llama_destination_id = accessor2->second;

#else
    try {
        scoped_lock<SpinLock> vlock(m_vmap_locks[e.m_source % m_num_vmap_locks]);
        llama_source_id = vmap_write_store_find(e.m_source); // throws std::out_of_range if ext_source_id is not present
    } catch( std::out_of_range& e ){
        return false; // the source does not exist
    }

    try {
        scoped_lock<SpinLock> vlock(m_vmap_locks[e.m_destination % m_num_vmap_locks]);
        llama_destination_id = vmap_write_store_find(e.m_destination); // throws std::out_of_range if ext_source_id is not present
    } catch( std::out_of_range& e ){
        return false; // the destination does not exist
    }
#endif

    // Again, this cannot be strictly correct, a race condition can still occur where a vertex is removed in the meanwhile we are
    // changing an edge
    if(!m_is_directed){ // for undirected graphs, ensure llama_source_id < llama_destination_id
        if(llama_source_id > llama_destination_id){
            std::swap(llama_source_id, llama_destination_id);
        }
    }

    return m_db->graph()->delete_edge_if_exists(llama_source_id, llama_destination_id);
}

void LLAMAClass::build(){
    scoped_lock<shared_mutex_t> xlock(m_lock_checkpoint);
#if defined(LL_PROFILE_UPDATES)
    auto t_start = chrono::steady_clock::now();
#endif

    COUT_DEBUG("build");

#if !defined(LLAMA_HASHMAP_WITH_TBB)
    // merge the changes to the vertex ids into m_vmap_read_only
    auto changeset_removed = m_vmap_removed.lock_table();
    for(auto it = begin(changeset_removed); it != end(changeset_removed); it++){
        m_vmap_read_only.erase(it->first);
    }
    changeset_removed.unlock();

    auto changeset_updated = m_vmap_updated.lock_table();
    for(auto it = begin(changeset_updated); it != end(changeset_updated); it++){
        m_vmap_read_only.insert(it->first, it->second);
    }
    changeset_updated.unlock();

    m_vmap_updated.clear();
    m_vmap_removed.clear();
#endif

    assert((static_cast<int64_t>(m_num_edges) + m_db->graph()->get_num_edges_diff() >= 0) && "Underflow");
    m_num_edges += m_db->graph()->get_num_edges_diff();

    // finally, create the new delta
    m_db->graph()->checkpoint();

#if defined(LL_PROFILE_UPDATES)
    common::compiler_barrier();
    g_llama_build_nanosecs += chrono::duration_cast<chrono::nanoseconds>(chrono::steady_clock::now() - t_start).count();
#endif
}


/*****************************************************************************
 *                                                                           *
 *  Start/stop the timer profiling                                            *
 *                                                                           *
 *****************************************************************************/
#if defined(LL_PROFILE_UPDATES)
void LLAMAClass::updates_start(){
    scoped_lock<shared_mutex_t> xlock(m_lock_checkpoint);
    g_llama_total_time.start();
}

void LLAMAClass::updates_stop(){
    scoped_lock<shared_mutex_t> xlock(m_lock_checkpoint);
    g_llama_total_time.stop();
}
#endif

/*****************************************************************************
 *                                                                           *
 *  Dump                                                                     *
 *                                                                           *
 *****************************************************************************/

void LLAMAClass::dump_ostream(std::ostream& out) const {
    shared_lock<shared_mutex_t> cplock(m_lock_checkpoint); // forbid any checkpoint now
    dump_impl(out, *(m_db->graph()));
}

void LLAMAClass::dump_ostream(std::ostream& out, int level) const {
    shared_lock<shared_mutex_t> cplock(m_lock_checkpoint); // forbid any checkpoint now

    // some sanity checks
    if(level < 0){
        throw std::invalid_argument("level < 0");
    } else if(level == (int) num_levels()){
        dump_impl(out, *(m_db->graph())); // write store
        return;
    } else if (level > (int) num_levels()){
        throw std::invalid_argument("level > num_levels()");
    }

    ll_mlcsr_ro_graph graph{ &(m_db->graph()->ro_graph()) , level };
    dump_impl(out, graph);
}

void LLAMAClass::dump_snapshot(ll_mlcsr_ro_graph& graph) const{
    dump_impl(std::cout, graph);
}

template<typename T>
void LLAMAClass::dump_impl(std::ostream& out, T& graph) const{
    auto names = reinterpret_cast<ll_mlcsr_node_property<uint64_t>*>(graph.get_node_property_64(g_llama_property_names));

    out << "[LLAMA] num vertices (global): " << num_vertices() << ", num edges (global): " << num_edges() << ", max nodes: " <<
            graph.max_nodes() << ", num levels (#snapshots): " << graph.num_levels() << endl;

    for(node_t node_id = 0; node_id < graph.max_nodes(); node_id++){
        // graph.node_exists is quite inconsistent. Depending whether the node is in the write store or in the read store, the behaviour is different
        // in the read store, the node exists if has at least one incoming or outgoing edge
        // in the write store, the node exists if it has just been inserted with #insert, regardless of the number of edges it has
        if(!(graph.node_exists(node_id))) continue;

        out << "[" << node_id << " -> " << names->get(node_id) << "] " << get_write_store_outdegree(node_id) << " outgoing edges: ";

        ll_edge_iterator iterator;
        graph.out_iter_begin(iterator, node_id);
        bool first = true;
        for (edge_t e = graph.out_iter_next(iterator); e != LL_NIL_EDGE; e = graph.out_iter_next(iterator)) {
            node_t n = LL_ITER_OUT_NEXT_NODE(graph, iterator, e);

            if(first) first = false; else out << ", ";
            if(e >= 0){
                out << "<" << names->get(n) << ", e_id: " << e << ", weight: " << get_out_edge_weight(graph, e) << ">";
            } else {
                out << names->get(n);
            }
        }

        if(!m_is_directed){
            graph.in_iter_begin_fast(iterator, node_id);
            for (edge_t e = graph.in_iter_next_fast(iterator); e != LL_NIL_EDGE; e = graph.in_iter_next_fast(iterator)) {
                node_t n = LL_ITER_OUT_NEXT_NODE(graph, iterator, e);

                if(first) first = false; else out << ", ";
                if(e >= 0){
                    out << "<" << names->get(n) << ", e_id: " << e << ", e_weight: " << get_in_edge_weight(graph, e) << ">";
                } else {
                    out << names->get(n);
                }
            }
        }

        out << "\n";
    }
}

// explicitly instantiate the method for both ll_mlcsr_ro_graph and ll_writable_graph
template
void LLAMAClass::dump_impl<ll_mlcsr_ro_graph>(std::ostream& out, ll_mlcsr_ro_graph& graph) const;
template
void LLAMAClass::dump_impl<ll_writable_graph>(std::ostream& out, ll_writable_graph& graph) const;


} // namespace
