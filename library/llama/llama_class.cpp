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

using namespace common;
using namespace std;

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG
extern mutex _log_mutex [[maybe_unused]];
#define COUT_DEBUG_FORCE(msg) { std::scoped_lock<std::mutex> lock{::_log_mutex}; std::cout << "[LLAMA::" << __FUNCTION__ << "] " << msg << std::endl; }
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
// sometimes I'm referring them as singular, sometimes as plural. Let's use a const for peace of mind
char const * const g_llama_property_names = "names";
char const * const g_llama_property_weights = "weights";

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

/*****************************************************************************
 *                                                                           *
 *  Init                                                                     *
 *                                                                           *
 *****************************************************************************/

namespace library {

LLAMAClass::LLAMAClass(bool is_directed) : m_is_directed(is_directed) {
    m_db = new ll_database();

    auto& csr = m_db->graph()->ro_graph();
    csr.create_uninitialized_node_property_64(g_llama_property_names, LL_T_INT64); // translate a logical id (e.g. node_id = 4) into an external node id (e.g. user_id = 21084718954)
    csr.create_uninitialized_edge_property_64(g_llama_property_weights, LL_T_DOUBLE);
}


LLAMAClass::~LLAMAClass(){
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
    return m_vmap_read_only.size();
}

uint64_t LLAMAClass::num_levels() const{
    return m_db->graph()->ro_graph().num_levels();
}

bool LLAMAClass::is_directed() const {
    return m_is_directed;
}

bool LLAMAClass::has_vertex(uint64_t vertex_id) const {
    try {
        vmap_write_store_find(vertex_id); // exception -> vertex not present
        return true;
    } catch (std::out_of_range& e){
        return false;
    }
}

double LLAMAClass::get_weight(uint64_t ext_source_id, uint64_t ext_destination_id) const {
    auto nan = numeric_limits<double>::quiet_NaN();

    uint64_t int_source_id, int_destination_id;
    try {
        scoped_lock<SpinLock> lock(m_lock_vertex_map);
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
}

void LLAMAClass::set_timeout(uint64_t seconds) {
    m_timeout = chrono::seconds{ seconds };
}

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
    COUT_DEBUG("vertex_id: " << vertex_id_ext);
    shared_lock<shared_mutex> cplock(m_lock_checkpoint); // forbid any checkpoint now

    m_lock_vertex_map.lock();
    if(!vmap_write_store_contains(vertex_id_ext)){
        m_lock_vertex_map.unlock();

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
        m_lock_vertex_map.unlock();
        return false;
    }
}

bool LLAMAClass::remove_vertex(uint64_t vertex_id_ext){
    COUT_DEBUG("vertex_id: " << vertex_id_ext);
    shared_lock<shared_mutex> cplock(m_lock_checkpoint); // forbid any checkpoint now

    bool is_removed = false;
    int64_t llama_vertex_id = -1;

    m_lock_vertex_map.lock();
    if(m_vmap_read_only.find(vertex_id_ext, llama_vertex_id)){
        m_vmap_removed.insert(vertex_id_ext, true);
        is_removed = true;
    }

    is_removed |= m_vmap_updated.erase_fn(vertex_id_ext, [&llama_vertex_id](int64_t& mapping){
        llama_vertex_id = mapping;
        return true;
    });


    m_lock_vertex_map.unlock();

    if(is_removed){
        assert(llama_vertex_id != -1);
        m_db->graph()->delete_node(llama_vertex_id);
    }

    return is_removed;
}

bool LLAMAClass::add_edge(graph::WeightedEdge e){
    COUT_DEBUG("edge: " << e);
    shared_lock<shared_mutex> cplock(m_lock_checkpoint); // forbid any checkpoint now

    node_t llama_source_id { -1 }, llama_destination_id { -1 };
    try {
        scoped_lock<SpinLock> lock(m_lock_vertex_map);
        llama_source_id = vmap_write_store_find(e.m_source); // throws std::out_of_range if ext_source_id is not present
        llama_destination_id = vmap_write_store_find(e.m_destination); // as above
    } catch( std::out_of_range& e ){
        return false; // either the source or the destination does not exist
    }

    if(!m_is_directed){ // for undirected graphs, ensure llama_source_id < llama_destination_id
        if(llama_source_id > llama_destination_id){
            std::swap(llama_source_id, llama_destination_id);
        }
    }

    edge_t edge_id;
    bool inserted = m_db->graph()->add_edge_if_not_exists(llama_source_id, llama_destination_id, &edge_id);

    // thread unsafe, this should really still be under the same latch of add_edge_if_not_exists
    if(inserted){
        m_db->graph()->get_edge_property_64(g_llama_property_weights)->add(edge_id, *reinterpret_cast<uint64_t*>(&(e.m_weight)));
    }

    return inserted;
}

bool LLAMAClass::remove_edge(graph::Edge e){
    COUT_DEBUG("edge: " << e);
    shared_lock<shared_mutex> cplock(m_lock_checkpoint); // forbid any checkpoint now

    node_t llama_source_id { -1 }, llama_destination_id { -1 };
    try {
        scoped_lock<SpinLock> lock(m_lock_vertex_map);
        llama_source_id = vmap_write_store_find(e.m_source); // throws std::out_of_range if ext_source_id is not present
        llama_destination_id = vmap_write_store_find(e.m_destination); // as above

        if(!m_is_directed){ // for undirected graphs, ensure llama_source_id < llama_destination_id
            if(llama_source_id > llama_destination_id){
                std::swap(llama_source_id, llama_destination_id);
            }
        }

//        llama_edge_id = m_db->graph()->find(llama_source_id, llama_destination_id);
    } catch(std::out_of_range& e){
        return false; // either source or destination do not exist
    }

    return m_db->graph()->delete_edge_if_exists(llama_source_id, llama_destination_id);

//    if(llama_edge_id != LL_NIL_EDGE){
////        m_db->graph()->delete_edge(llama_source_id, llama_edge_id);
//        m_db-
//        return true;
//    } else {
//        return false;
//    }
}

void LLAMAClass::build(){
    COUT_DEBUG("build");
    scoped_lock<shared_mutex> xlock(m_lock_checkpoint);

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

    assert((static_cast<int64_t>(m_num_edges) + m_db->graph()->get_num_edges_diff() >= 0) && "Underflow");
    m_num_edges += m_db->graph()->get_num_edges_diff();

    // finally, create the new delta
    m_db->graph()->checkpoint();
    m_db->graph()->callback_ro_changed();
}


/*****************************************************************************
 *                                                                           *
 *  Dump                                                                     *
 *                                                                           *
 *****************************************************************************/

void LLAMAClass::dump_ostream(std::ostream& out) const {
    shared_lock<shared_mutex> cplock(m_lock_checkpoint); // forbid any checkpoint now
    dump_impl(out, *(m_db->graph()));
}

void LLAMAClass::dump_ostream(std::ostream& out, int level) const {
    shared_lock<shared_mutex> cplock(m_lock_checkpoint); // forbid any checkpoint now

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
