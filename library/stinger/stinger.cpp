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

// stdlib & system includes
#include <cassert>
#include <cstdlib> // strtoull
#include <mutex>
#include <sstream>
#include <string>

// stinger support
#include "stinger_core/stinger.h"
#include "stinger_utils/metisish_support.h"

// internal includes
#include "reader/format.hpp" // graph format for the method #load

using namespace std;

namespace library {

// Macros
#define STINGER reinterpret_cast<struct stinger*>(m_stinger_graph)
#undef CURRENT_ERROR_TYPE
#define CURRENT_ERROR_TYPE ::library::StingerError

/******************************************************************************
 *                                                                            *
 *  Debug                                                                     *
 *                                                                            *
 *****************************************************************************/
static mutex _stinger_debug_mutex;
//#define DEBUG
#define COUT_DEBUG_FORCE(msg) { scoped_lock<mutex> lock(_stinger_debug_mutex); std::cout << "[Stinger::" << __FUNCTION__ << "] [" << this_thread::get_id() << "] " << msg << std::endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif

/******************************************************************************
 *                                                                            *
 *  Initialisation                                                            *
 *                                                                            *
 *****************************************************************************/

Stinger::Stinger(){
    m_stinger_graph = stinger_new();
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
// Covert a vertex id from the external domain to the stinger internal id
// @return a value >= 0 than if it's the index in the stinger adjacency list, or -1 in case of error
static int64_t convert_external2stinger(struct stinger* graph, uint64_t external_vertex_id){
    char buffer[64];
    sprintf(buffer, "%" PRIu64, external_vertex_id);
    return stinger_mapping_lookup(graph, buffer, sizeof(buffer));
}

/******************************************************************************
 *                                                                            *
 *  Properties                                                                *
 *                                                                            *
 *****************************************************************************/
uint64_t Stinger::num_edges() const {
    return static_cast<uint64_t>(stinger_total_edges(STINGER));
}

uint64_t Stinger::num_vertices() const {
    return stinger_num_active_vertices(STINGER);
}

bool Stinger::has_vertex(uint64_t vertex_id) const {
    char buffer[64];
    sprintf(buffer, "%" PRIu64, vertex_id);
    int64_t stinger_vertex_id = stinger_mapping_lookup(STINGER, buffer, sizeof(buffer)); // stinger_names is reported to be thread_safe
    return stinger_vertex_id >= 0; // -1 => the vertex_id has not been set in the dictionary
}

bool Stinger::has_edge(uint64_t ext_source_id, uint64_t ext_destination_id) const {
    char buffer[64];
    sprintf(buffer, "%" PRIu64, ext_source_id);
    int64_t stg_source_id = stinger_mapping_lookup(STINGER, buffer, sizeof(buffer)); // stinger_names is reported to be thread_safe
    if(stg_source_id <= 0) return false; // the source node does not exist
    sprintf(buffer, "%" PRIu64, ext_destination_id);
    int64_t stg_destination_id = stinger_mapping_lookup(STINGER, buffer, sizeof(buffer)); // stinger_names is reported to be thread_safe
    if(stg_destination_id <= 0) return false; // the destination node does not exist

    bool exists = false;
    STINGER_FORALL_OUT_EDGES_OF_VTX_BEGIN(STINGER, stg_source_id) {
        if(STINGER_EDGE_DEST == stg_destination_id){
            exists = true;
        }
    } STINGER_FORALL_OUT_EDGES_OF_VTX_END();

    return exists;
}


/******************************************************************************
 *                                                                            *
 *  Update                                                                    *
 *                                                                            *
 *****************************************************************************/
bool Stinger::add_vertex(uint64_t vertex_id){
    // Stinger has an internal dictionary where the vertex_id from the external are mapped
    // into internal indices for the internal adjacency list. Here, we simply register
    // the mapping to the dictionary.

    char buffer[65];
    sprintf (buffer, "%" PRIu64, vertex_id);
    int64_t out_stg_vertex_id { -1 }; // ignored
    // @return 0 on existing mapping found. 1 on new mapping created. -1 on failure (STINGER is full).
    int rc = stinger_mapping_create(STINGER, buffer, sizeof(buffer), &out_stg_vertex_id);

    switch(rc){
    case -1:
        ERROR("Cannot insert the given vertex_id: " << vertex_id << ". Stinger is full");
    case 0:
        return false; // the vertex_id already exists
    case 1:
        return true; // mapping created
    default:
        ERROR("Invalid return code: " << rc);
    }
}

bool Stinger::delete_vertex(uint64_t ext_vertex_id){
    int64_t vertex_id = convert_external2stinger(STINGER, ext_vertex_id);
    if(vertex_id < 0) return false; // the vertex does not exist

    auto rc = stinger_remove_vertex(STINGER, ext_vertex_id);
    // the implementation already removes the vertex_id from the  internal naming dictionary
    // @return 0 if deletion is successful, -1 if it fails
    COUT_DEBUG("vertex_id: " << ext_vertex_id << " (internal: " << vertex_id << "), rc: " << rc);

    return rc == 0;
}

bool Stinger::add_edge(graph::WeightedEdge e){
    // get the indices in the adjacency lists
    int64_t src = convert_external2stinger(STINGER, e.source());
    if(src < 0) return false; // the source does not exist
    int64_t dst = convert_external2stinger(STINGER, e.destination());
    if(dst < 0) return false; // the destination does not exist

    int rc = stinger_insert_edge (STINGER, /* type, ignore */ 0 , src, dst, static_cast<int64_t>(e.weight()), /* timestamp, ignore */ 0);
    // @return 1 if edge is inserted successfully for the first time, 0 if edge is already found and updated, -1 if error.
    COUT_DEBUG("source: " << e.source() << " (internal: " << src << "), destination: " << e.destination() << " (internal: " << dst << "), weight: " << e.weight() << ", rc: " << rc);

    return (rc == 0 || rc == 1);
}

bool Stinger::delete_edge(graph::Edge e){
    // get the indices in the adjacency lists
    int64_t src = convert_external2stinger(STINGER, e.source());
    if(src < 0) return false; // the source does not exist
    int64_t dst = convert_external2stinger(STINGER, e.destination());
    if(dst < 0) return false; // the destination does not exist

    int rc = stinger_remove_edge(STINGER, /* type, ignore */ 0, src, dst);
    // @return 1 on success, 0 if the edge is not found.
    COUT_DEBUG("source: " << e.source() << " (internal: " << src << "), destination: " << e.destination() << " (internal: " << dst << "), rc: " << rc);

    return (rc == 1);
}

/******************************************************************************
 *                                                                            *
 *  Load                                                                      *
 *                                                                            *
 *****************************************************************************/

void Stinger::load(const std::string& path){
    using namespace reader;
    auto graph_format = get_graph_format(path);
    switch(graph_format){
    case Format::METIS: {
        // it always returns 0, in case of error it just terminates the program on its own
        load_metisish_graph(STINGER, (char*) path.c_str()); // safe to cast away the constant qualifier
    } break;
    default:
        ERROR("Graph format not supported: " << path);
    }
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

static void dump(struct stinger* graph){
    cout << "Vertices: " << stinger_num_active_vertices(graph) << ", edges: " << stinger_total_edges(graph) << ", size: " << stinger_graph_size(graph) << " bytes" << "\n";
    for(int64_t vertex_id = 0, vertex_max = stinger_max_active_vertex(graph); vertex_id <= vertex_max; vertex_id++){
        if(stinger_degree_get(graph, vertex_id) == 0) continue; // gap, empty slot not used

        char* vertex_id_name = nullptr; uint64_t vertex_id_name_length = 0;
        int rc = stinger_mapping_physid_get(graph, vertex_id, &vertex_id_name, &vertex_id_name_length);
        if(rc == 0){
            cout << "[" << vertex_id_name << ", internal ID: " << vertex_id << "] ";
        } else {
            cout << "[" << vertex_id << "] ";
        }
        free(vertex_id_name); vertex_id_name = nullptr;

        cout << "type: " << dump_vertex_type(graph, vertex_id) << ", weight: " << stinger_vweight_get(graph, vertex_id) << ", " <<
                "degree in/out/total: " << stinger_indegree_get(graph, vertex_id) << "/" << stinger_outdegree_get(graph, vertex_id) << "/" << stinger_degree_get(graph, vertex_id);
        if(stinger_outdegree_get(graph, vertex_id) == 0){
            cout << "\n";
        } else {
            cout << ", outgoing edges: \n";
            STINGER_FORALL_OUT_EDGES_OF_VTX_BEGIN(graph, vertex_id) {
                cout << "  ";
                char* vertex_id_name = nullptr; uint64_t vertex_id_name_length = 0;
                int rc = stinger_mapping_physid_get(graph, STINGER_EDGE_DEST, &vertex_id_name, &vertex_id_name_length);
                if(rc == 0){
                    cout << vertex_id_name << ", internal ID: " << STINGER_EDGE_DEST;
                } else {
                    cout << STINGER_EDGE_DEST;
                }
                free(vertex_id_name); vertex_id_name = nullptr;
                cout << ", ";

                cout << "type: " << dump_edge_type(graph, STINGER_EDGE_TYPE) << ", " <<
                        "weight: " << STINGER_EDGE_WEIGHT << "\n";
            } STINGER_FORALL_OUT_EDGES_OF_VTX_END();
        }
    }
}

void Stinger::dump() const {
    ::library::dump(STINGER); // trampoline
}

} // namespace
