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

#include <cassert>
#include <mutex>

#include "library/stinger/stinger_error.hpp"
#include "stinger_core/stinger.h"

using namespace common;
using namespace libcuckoo;
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
#define COUT_DEBUG_FORCE(msg) { std::scoped_lock<std::mutex> lock{::_log_mutex}; std::cout << "[Stinger::" << __FUNCTION__ << "] " << msg << std::endl; }
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

/*****************************************************************************
 *                                                                           *
 *  Initialisation                                                           *
 *                                                                           *
 *****************************************************************************/
namespace gfe::library {

Stinger::Stinger(bool directed) : m_directed(directed){
    struct stinger_config_t config;
    memset(&config, 0, sizeof(config)); // init
    config.nv = 1ll<<32; // max number of vertices, 4G
    config.nebs = 0; // max number of edge blocks, 0=auto
    config.netypes = 1; // number of edge types, we are not going to use this feature, 1
#if defined(STINGER_USE_INTERNAL_MAPPING)
    config.nvtypes = 2; // number of vertex types, 0=vertex active, 1=vertex removed, 2
#else
    config.nvtypes = 1; // number of vertex types, we are not going to use this feature, 1
#endif
    config.no_map_none_etype = config.no_map_none_vtype = 1; // I think this is whether we want to attach names (as strings) to the vertex/edge types, 1 = we don't use this feature
    config.no_resize = 0; // when building stinger for the first time, allow to resize the internal structs if they do not fit into memory, 0 = allow resize

    m_stinger_graph = stinger_new_full(&config);
    assert(m_stinger_graph != nullptr && "Stinger allocation");
    if(m_stinger_graph == nullptr) ERROR("Cannot allocate stinger graph");
}

Stinger::~Stinger(){
    stinger_free_all(STINGER); m_stinger_graph = nullptr;
}

/******************************************************************************
 *                                                                            *
 *  Conversion functions                                                      *
 *                                                                            *
 *****************************************************************************/
#if defined(STINGER_USE_INTERNAL_MAPPING)
// Covert a vertex id from the external domain to the stinger internal id
// @return a value >= 0 than if it's the index in the stinger adjacency list, or -1 in case of error
static int64_t convert_external2stinger(struct stinger* graph, uint64_t external_vertex_id){
    char buffer[64];
    sprintf(buffer, "%" PRIu64, external_vertex_id);
    return stinger_mapping_lookup(graph, buffer, sizeof(buffer));
}

int64_t Stinger::get_internal_id(uint64_t vertex_id) const{
    return convert_external2stinger(STINGER, vertex_id);
}

uint64_t Stinger::get_external_id(int64_t internal_vertex_id) const{
    uint64_t result = numeric_limits<uint64_t>::max();

    char* vertex_id_name = nullptr; uint64_t vertex_id_name_length = 0;
    int rc = stinger_mapping_physid_get(STINGER, internal_vertex_id, &vertex_id_name, &vertex_id_name_length);
    if( rc == 0 && vertex_id_name_length > 0 ){ // mapping found
        result = stoull(vertex_id_name);
    }

    free(vertex_id_name); vertex_id_name = nullptr;

    return result;
}

uint64_t Stinger::get_max_num_mappings() const {
    return stinger_mapping_nv(STINGER);
}

template<typename T>
void Stinger::to_external_ids(const T* __restrict internal_ids, size_t internal_ids_sz, cuckoohash_map<uint64_t, T>* external_ids){
    ASSERT(external_ids != nullptr);

    #pragma omp parallel for
    for(uint64_t internal_id = 0; internal_id < internal_ids_sz; internal_id++){
        if(stinger_vtype_get(STINGER, internal_id) == 0){ // if = 1, the node is marked for deletion
            char* vertex_id_name = nullptr; uint64_t vertex_id_name_length = 0;
            int rc = stinger_mapping_physid_get(STINGER, internal_id, &vertex_id_name, &vertex_id_name_length);
            if( rc == 0 ){ // mapping found
                uint64_t external_id = stoull(vertex_id_name);
                COUT_DEBUG("external_id: " << external_id << ", internal_id: " << internal_id);
                external_ids->insert(external_id, internal_ids[internal_id]);
            }
            free(vertex_id_name); vertex_id_name = nullptr;
        }
    }
}

int64_t Stinger::get_or_create_vertex_id(uint64_t external_vertex_id) {
    char buffer[64];
    sprintf(buffer, "%" PRIu64, external_vertex_id);

    int64_t internal_vertex_id { -1 }; // output
    // @return 0 on existing mapping found. 1 on new mapping created. -1 on failure (STINGER is full).
    int rc = stinger_mapping_create(STINGER, buffer, sizeof(buffer), &internal_vertex_id);

    switch(rc){
    case -1:
        ERROR("Cannot insert the given vertex_id: " << external_vertex_id << ". Stinger is full");
    case 0: { // the mapping already exists
        return internal_vertex_id;
    } break;
    case 1: { // mapping created
        m_num_vertices++;
        return internal_vertex_id; // mapping created
    }
    default:
        ERROR("Invalid return code: " << rc);
    }
}

#else
int64_t Stinger::get_internal_id(uint64_t external_vertex_id) const {
    try {
        return m_vertex_mappings_e2i.find(external_vertex_id);
    } catch(std::out_of_range& e){ // not found
        return -1;
    }
}

uint64_t Stinger::get_external_id(int64_t internal_vertex_id) const{
    try {
        return m_vertex_mappings_i2e.find(internal_vertex_id);
    } catch(std::out_of_range& e){ // not found
        return numeric_limits<uint64_t>::max();
    }
}

uint64_t Stinger::get_max_num_mappings() const {
    scoped_lock<SpinLock> lock(m_spin_lock);
    return m_next_vertex_id;
}

template<typename T>
void Stinger::to_external_ids(const T* __restrict internal_ids, size_t internal_ids_sz, cuckoohash_map<uint64_t, T>* external_ids){
    ASSERT(external_ids != nullptr);

    #pragma omp parallel for
    for(uint64_t internal_id = 0; internal_id < internal_ids_sz; internal_id++){
        uint64_t external_id = get_external_id(internal_id);
//        COUT_DEBUG("internal_id: " << internal_id << " -> " << external_id);
        if(external_id != numeric_limits<uint64_t>::max()){
            external_ids->insert(external_id, internal_ids[internal_id]);
        }
    }
}
#endif

template<typename T>
void Stinger::to_external_ids(const vector<T>& internal_ids, cuckoohash_map<uint64_t, T>* external_ids){
    to_external_ids(internal_ids.data(), internal_ids.size(), external_ids);
}

template void Stinger::to_external_ids<double>(const vector<double>& internal_ids, cuckoohash_map<uint64_t, double>* external_ids);
template void Stinger::to_external_ids<double>(const double* __restrict internal_ids, size_t internal_ids_sz, cuckoohash_map<uint64_t, double>* external_ids);
template void Stinger::to_external_ids<int64_t>(const vector<int64_t>& internal_ids, cuckoohash_map<uint64_t, int64_t>* external_ids);
template void Stinger::to_external_ids<int64_t>(const int64_t* __restrict internal_ids, size_t internal_ids_sz, cuckoohash_map<uint64_t, int64_t>* external_ids);

/******************************************************************************
 *                                                                            *
 *  Properties                                                                *
 *                                                                            *
 *****************************************************************************/
uint64_t Stinger::num_edges() const {
    uint64_t num_directed_edges = static_cast<uint64_t>(stinger_total_edges(STINGER));
    if(!m_directed) num_directed_edges /= 2; /* divided by two, because we always store the edge for both source/destination */
    return num_directed_edges;
}

uint64_t Stinger::num_vertices() const {
    return m_num_vertices;
}

bool Stinger::has_vertex(uint64_t external_vertex_id) const {
    int64_t internal_vertex_id = get_internal_id(external_vertex_id);
    if(internal_vertex_id < 0) return false; // the mapping does not exist

#if defined(STINGER_USE_INTERNAL_MAPPING)
    return stinger_vtype_get(STINGER, internal_vertex_id) == 0;
#else
    return true;
#endif
}

double Stinger::get_weight(uint64_t ext_source_id, uint64_t ext_destination_id) const {
    COUT_DEBUG("external: " << ext_source_id << " -> " << ext_destination_id);
    constexpr double NaN { numeric_limits<double>::signaling_NaN() };
    int64_t int_source_id = get_internal_id(ext_source_id);
    if(int_source_id < 0) return NaN; // the mapping does not exist
    int64_t int_destination_id = get_internal_id(ext_destination_id);
    if(int_destination_id < 0) return NaN; // the mapping does not exist

    STINGER_FORALL_OUT_EDGES_OF_VTX_BEGIN(STINGER, int_source_id) {
        if(STINGER_EDGE_DEST == int_destination_id){
            return INT2DBL(STINGER_EDGE_WEIGHT);
        }
    } STINGER_FORALL_OUT_EDGES_OF_VTX_END();

    return NaN;
}

bool Stinger::is_directed() const {
    return m_directed;
}

void Stinger::set_timeout(uint64_t seconds) {
    m_timeout = seconds;
}

void* Stinger::handle(){
    return m_stinger_graph;
}

/******************************************************************************
 *                                                                            *
 *  Updates                                                                   *
 *                                                                            *
 *****************************************************************************/
#if defined(STINGER_USE_INTERNAL_MAPPING)
bool Stinger::add_vertex(uint64_t vertex_id){
    COUT_DEBUG("vertex_id: " << vertex_id);

    // Stinger has an internal dictionary where the vertex_id from the external are mapped
    // into internal indices for the internal adjacency list. Here, we simply register
    // the mapping to the dictionary.

    char buffer[65];
    sprintf (buffer, "%" PRIu64, vertex_id);
    int64_t internal_vertex_id { -1 }; // output
    // @return 0 on existing mapping found. 1 on new mapping created. -1 on failure (STINGER is full).
    int rc = stinger_mapping_create(STINGER, buffer, sizeof(buffer), &internal_vertex_id);

    switch(rc){
    case -1:
        ERROR("Cannot insert the given vertex_id: " << vertex_id << ". Stinger is full");
    case 0: // the mapping already exists
        return false;
//    case 0: {
//        // because we cannot delete an existing mapping, here we use the convention of setting a vertex type to 1
//        // if the vertex is supposed to be deleted
//        scoped_lock<SpinLock> lock(m_spin_lock);
//
//        vtype_t type = stinger_vtype_get(STINGER, internal_vertex_id);
//        switch(type){
//        case 0:// the mapping already existed and the vertex is active (type = 0)
//            return false;
//        case 1: // the mapping already existed but the vertex was deleted (type = 1)
//            stinger_vtype_set(STINGER, internal_vertex_id, 0); // reset to active
//            m_num_vertices ++;
//            return true;
//        default:
//            ERROR("Invalid type: " << type);
//        }
//    } break;
    case 1: {
        m_num_vertices++;
        return true; // mapping created
    }
    default:
        ERROR("Invalid return code: " << rc);
    }
}

bool Stinger::remove_vertex(uint64_t external_vertex_id){
    // remove all edges (directed and undirected) for the internal_vertex_id.
    // It should also remove the mapping, but the current Stinger implementation does not support this, returning -1.
    //auto rc = stinger_remove_vertex(STINGER, internal_vertex_id);

    assert(0 && "Operation not supported");
    return false;

//    int64_t internal_vertex_id = get_internal_id(external_vertex_id);
//    if(internal_vertex_id < 0) return false; // the vertex does not even have a mapping
//
//    // remove all edges (directed and undirected) for the internal_vertex_id.
//    // It should also remove the mapping, but the current Stinger implementation does not support this, returning -1.
//    // As workaround, we set the vertex type to 1, to signal the vertex as inactive
//    auto rc = stinger_remove_vertex(STINGER, internal_vertex_id);
//
//    scoped_lock<SpinLock> lock(m_spin_lock);
//    if(rc == -1){
//        if(stinger_vtype_get(STINGER, internal_vertex_id) == 1){ // already previously marked as inactive
//            return false;
//        } else {
//            stinger_vtype_set(STINGER, internal_vertex_id, 1);
//        }
//    }
//    assert(m_num_vertices > 0);
//    m_num_vertices--;

    return true;
}
#else
bool Stinger::add_vertex(uint64_t vertex_id){
    COUT_DEBUG("vertex_id: " << vertex_id);

    scoped_lock<SpinLock> lock(m_spin_lock);
    if(m_vertex_mappings_e2i.contains(vertex_id)){
        return false;
    } else {
        int64_t mapping;
        if(!m_reuse_vertices.empty()){
            mapping = m_reuse_vertices.back(); m_reuse_vertices.pop_back();
        } else {
            mapping = m_next_vertex_id++;
        }

        m_vertex_mappings_i2e.insert(mapping, vertex_id);
        m_vertex_mappings_e2i.insert(vertex_id, mapping);
        m_num_vertices++;

        return true;
    }
}

bool Stinger::remove_vertex(uint64_t external_vertex_id){
    int64_t internal_vertex_id = get_internal_id(external_vertex_id);
    if(internal_vertex_id < 0) return false; // the vertex does not even have a mapping

    scoped_lock<SpinLock> lock(m_spin_lock);

    m_vertex_mappings_e2i.erase(external_vertex_id);
    m_vertex_mappings_i2e.erase(internal_vertex_id);

    assert(m_num_vertices > 0);
    m_num_vertices--;
    m_reuse_vertices.push_back(internal_vertex_id);

    return true;
}
#endif

bool Stinger::add_edge(graph::WeightedEdge e){
    COUT_DEBUG("edge: " << e);

    // get the indices in the adjacency lists
    int64_t src = get_internal_id(e.source());
    if(src < 0) return false; // the source does not exist
    int64_t dst = get_internal_id(e.destination());
    if(dst < 0) return false; // the destination does not exist
    int64_t weight = DBL2INT(e.m_weight);

    int rc = 0;

    if(m_directed){ // directed graph
        rc = stinger_insert_edge (STINGER, /* type, ignore */ 0 , src, dst, weight, /* timestamp, ignore */ 0);
    } else { // undirected graph
        rc = stinger_insert_edge_pair(STINGER, /* type, ignore */ 0, src, dst, weight, /* timestamp, ignore */ 0);
    }

    return rc >= 0;
}

bool Stinger::add_edge_v2(gfe::graph::WeightedEdge e){
    COUT_DEBUG("edge: " << e);

    // get the indices in the adjacency lists
    int64_t src = get_or_create_vertex_id(e.source());
    int64_t dst = get_or_create_vertex_id(e.destination());
    int64_t weight = DBL2INT(e.m_weight);

    int rc = 0;
    if(m_directed){ // directed graph
        rc = stinger_insert_edge (STINGER, /* type, ignore */ 0 , src, dst, weight, /* timestamp, ignore */ 0);
    } else { // undirected graph
        rc = stinger_insert_edge_pair(STINGER, /* type, ignore */ 0, src, dst, weight, /* timestamp, ignore */ 0);
    }

    return rc >= 0;
}

bool Stinger::remove_edge(graph::Edge e){
    // get the indices in the adjacency lists
    int64_t src = get_internal_id(e.source());
    if(src < 0) return false; // the source does not exist
    int64_t dst = get_internal_id(e.destination());
    if(dst < 0) return false; // the destination does not exist

    int rc = 0;
    if(m_directed){ // directed graph
        // @return 1 on success, 0 if the edge is not found.
        rc = stinger_remove_edge(STINGER, /* type, ignore */ 0, src, dst);

    } else { // undirected graph
        // @return < 0 in case of error, 0 edges not found, > 0 edges removed
        rc = stinger_remove_edge_pair(STINGER, /* type, ignore */ 0, src, dst);
    }

    return rc > 0;
}

/******************************************************************************
 *                                                                            *
 *  Dump                                                                      *
 *                                                                            *
 *****************************************************************************/
static string dump_vertex_type(struct stinger* graph, uint64_t vertex_id){
    stringstream ss;
    vtype_t vtype = stinger_vtype_get(graph, vertex_id);
    char* vname = stinger_vtype_names_lookup_name(graph, vtype);
    if(vname != nullptr && vname[0] != '\0'){
        ss << vname << " (" << vtype << ")";
    } else {
        ss << vtype;
    }
    return ss.str();
}

static string dump_edge_type(struct stinger* graph, int64_t edge_type){
    stringstream ss;
    char* edge_type_name = stinger_etype_names_lookup_name(graph, edge_type);
    if(edge_type_name != nullptr && edge_type_name[0] != '\0'){
        ss << edge_type_name << " (" << edge_type << ")";
    } else {
        ss << edge_type;
    }
    return ss.str();
}

void Stinger::dump_ostream(ostream& out) const {
    struct stinger* graph = STINGER;

    out << "[STINGER] Vertices: " << num_vertices() << ", edges: " << num_edges() << ", directed: " << std::boolalpha << is_directed() << ", size: " << stinger_graph_size(graph) << " bytes" << "\n";
    for(int64_t internal_vertex_id = 0, vertex_max = stinger_max_active_vertex(graph); internal_vertex_id <= vertex_max; internal_vertex_id++){
        if(stinger_vtype_get(graph, internal_vertex_id) == 1) continue; // vertex flagged as inactive

        uint64_t external_vertex_id = get_external_id(internal_vertex_id);
        if(external_vertex_id != numeric_limits<uint64_t>::max()){
            out << "[" << external_vertex_id << " (internal ID: " << internal_vertex_id << ")] ";
        } else {
            out << "[" << internal_vertex_id << "] ";
        }

        out << "type: " << dump_vertex_type(graph, internal_vertex_id) << ", weight: " << stinger_vweight_get(graph, internal_vertex_id) << ", " <<
                "degree in/out/total: " << stinger_indegree_get(graph, internal_vertex_id) << "/" << stinger_outdegree_get(graph, internal_vertex_id) << "/" << stinger_degree_get(graph, internal_vertex_id);
        if(stinger_outdegree_get(graph, internal_vertex_id) == 0){
            out << "\n";
        } else {
            out << ", outgoing edges: \n";
            STINGER_FORALL_OUT_EDGES_OF_VTX_BEGIN(graph, internal_vertex_id) {
                out << "  ";

                uint64_t vertex_id_name = get_external_id(STINGER_EDGE_DEST);
                if(vertex_id_name != numeric_limits<uint64_t>::max()){
                    out << vertex_id_name << ", internal ID: " << STINGER_EDGE_DEST;
                } else {
                    out << STINGER_EDGE_DEST;
                }
                out << ", ";

                out << "type: " << dump_edge_type(graph, STINGER_EDGE_TYPE) << ", " <<
                        "weight: " << INT2DBL(STINGER_EDGE_WEIGHT) << "\n";
            } STINGER_FORALL_OUT_EDGES_OF_VTX_END();
        }
    }
    std::flush(out);
}

} // namespace
