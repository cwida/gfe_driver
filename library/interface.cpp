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

#include "interface.hpp"

#include <cassert>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <mutex>

#include "common/error.hpp"
#include "common/quantity.hpp"
#include "common/system.hpp"
#include "common/timer.hpp"

#include "baseline/adjacency_list.hpp"
#include "baseline/dummy.hpp"
#if defined(HAVE_LLAMA)
#include "llama/llama_class.hpp"
#include "llama/llama_ref.hpp"
#include "llama-dv/llama-dv.hpp"
#endif
#include "reader/reader.hpp"
#if defined(HAVE_STINGER)
#include "stinger/stinger.hpp"
#include "stinger-dv/stinger-dv.hpp" // dense domain of vertices
#endif
#if defined(HAVE_GRAPHONE)
#include "graphone/graphone.hpp"
#endif
#if defined(HAVE_TESEO)
#include "teseo/teseo_driver.hpp"
#endif

using namespace std;

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG
namespace gfe { extern mutex _log_mutex [[maybe_unused]]; }
#define COUT_DEBUG_FORCE(msg) { std::scoped_lock<std::mutex> lock{::gfe::_log_mutex}; std::cout << "[Interface::" << __FUNCTION__ << "] " << msg << std::endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif


namespace gfe::library {
/*****************************************************************************
 *                                                                           *
 *  Factory                                                                  *
 *                                                                           *
 *****************************************************************************/
ImplementationManifest::ImplementationManifest(const string& name, const string& description, unique_ptr<Interface> (*factory)(bool)) :
    m_name(name), m_description(description), m_factory(factory){ }

std::unique_ptr<Interface> generate_baseline_adjlist(bool directed_graph){ // directed or undirected graph
    return unique_ptr<Interface>{ new AdjacencyList(directed_graph) };
}

std::unique_ptr<Interface> generate_dummy(bool directed_graph){
    return unique_ptr<Interface>{ new Dummy(directed_graph) };
}

#if defined(HAVE_LLAMA)
std::unique_ptr<Interface> generate_llama(bool directed_graph){
    return unique_ptr<Interface>{ new LLAMAClass(directed_graph) };
}
std::unique_ptr<Interface> generate_llama_dv(bool directed_graph){
    return unique_ptr<Interface>{ new LLAMA_DV(directed_graph) };
}
std::unique_ptr<Interface> generate_llama_dv_nobw(bool directed_graph){
    return unique_ptr<Interface>{ new LLAMA_DV(directed_graph, /* blind writes */ false) };
}
std::unique_ptr<Interface> generate_llama_ref(bool directed_graph){
    return unique_ptr<Interface>{ new LLAMARef(directed_graph) };
}
#endif

#if defined(HAVE_STINGER)
std::unique_ptr<Interface> generate_stinger(bool directed_graph){
    return unique_ptr<Interface>{ new Stinger(directed_graph) };
}
std::unique_ptr<Interface> generate_stinger_dv(bool directed_graph){
    return unique_ptr<Interface>{ new StingerDV(directed_graph) };
}
std::unique_ptr<Interface> generate_stinger_ref(bool directed_graph){
    return unique_ptr<Interface>{ new StingerRef(directed_graph) };
}
#endif

#if defined(HAVE_GRAPHONE)
// Heuristics to set the max number of vertices that the implementation can sustain.
// We wish to create up to 4G vertices, but unfortunately this limit is too high for usage in a workstation machines,
// so we decrease the value depending on the amount of available memory
static uint64_t graphone_max_num_vertices(){
    uint64_t ram = common::get_total_ram(); // in bytes
    // Ensure that the vertex array does not use more than 20% of the available RAM
    uint64_t num_vertices = std::min<uint64_t>(4ull<<30, ram/(5 * 8)); // 8 bytes per vertex, 1/5 of the total ram

    { // critical section, for the output
        std::scoped_lock<std::mutex> lock{_log_mutex};
        cout << "GraphOne: capacity of the vertex array implicitly set to: " << common::ComputerQuantity(num_vertices) << " vertices " << endl;
    }

    return num_vertices;
}

std::unique_ptr<Interface> generate_graphone_cons_sp(bool directed_graph){
    uint64_t N = graphone_max_num_vertices();
    return unique_ptr<Interface>{ new GraphOne(directed_graph, /* vtx dict ? */ true, /* blind writes ? */ false, /* ignore build ? */ false, /* ref impl ? */ false, N) };
}
std::unique_ptr<Interface> generate_graphone_cons_dv(bool directed_graph){
    uint64_t N = graphone_max_num_vertices();
    return unique_ptr<Interface>{ new GraphOne(directed_graph, /* vtx dict ? */ false, /* blind writes ? */ false, /* ignore build ? */ false, /* ref impl ? */ false, N) };
}
std::unique_ptr<Interface> generate_graphone_bw_sp(bool directed_graph){
    uint64_t N = graphone_max_num_vertices();
    return unique_ptr<Interface>{ new GraphOne(directed_graph, /* vtx dict ? */ true, /* blind writes ? */ true, /* ignore build ? */ false, /* ref impl ? */ false,  N) };
}
std::unique_ptr<Interface> generate_graphone_bw_dv(bool directed_graph){
    uint64_t N = graphone_max_num_vertices();
    return unique_ptr<Interface>{ new GraphOne(directed_graph, /* vtx dict ? */ false, /* blind writes ? */ true, /* ignore build ? */ false, /* ref impl ? */ false, N) };
}
std::unique_ptr<Interface> generate_graphone_bw_sp_ignore_build(bool directed_graph){
    uint64_t N = graphone_max_num_vertices();
    return unique_ptr<Interface>{ new GraphOne(directed_graph, /* vtx dict ? */ true, /* blind writes ? */ true, /* ignore build ? */ true, /* ref impl ? */ false, N) };
}
std::unique_ptr<Interface> generate_graphone_bw_dv_ignore_build(bool directed_graph){
    uint64_t N = graphone_max_num_vertices();
    return unique_ptr<Interface>{ new GraphOne(directed_graph, /* vtx dict ? */ false, /* blind writes ? */ true, /* ignore build ? */ true, /* ref impl ? */ false, N) };
}
std::unique_ptr<Interface> generate_graphone_ref(bool directed_graph){
    uint64_t N = graphone_max_num_vertices();
    return unique_ptr<Interface>{ new GraphOne(directed_graph, /* vtx dict ? */ true, /* blind writes ? */ true, /* ignore build ? */ false, /* ref impl ? */ true, N) };
}
std::unique_ptr<Interface> generate_graphone_ref_ignore_build(bool directed_graph){
    uint64_t N = graphone_max_num_vertices();
    return unique_ptr<Interface>{ new GraphOne(directed_graph, /* vtx dict ? */ true, /* blind writes ? */ true, /* ignore build ? */ true, /* ref impl ? */ true, N) };
}
#endif

#if defined(HAVE_TESEO)
std::unique_ptr<Interface> generate_teseo(bool directed_graph){
    return unique_ptr<Interface>{ new TeseoDriver(directed_graph) };
}
#endif

vector<ImplementationManifest> implementations() {
    vector<ImplementationManifest> result;

    result.emplace_back("baseline", "Sequential baseline, based on adjacency list", &generate_baseline_adjlist);
    result.emplace_back("dummy", "Dummy implementation of the interface, all operations are nop", &generate_dummy);

#if defined(HAVE_LLAMA)
    // v2 25/11/2019: better scalability for the llama dictionary
    result.emplace_back("llama3", "LLAMA library", &generate_llama);
    result.emplace_back("llama-dv", "LLAMA with dense vertices", &generate_llama_dv);
    result.emplace_back("llama-dv-nobw", "LLAMA with dense vertices, no blind writes", &generate_llama_dv_nobw);
    result.emplace_back("llama3-ref", "LLAMA with the GAPBS ref impl.", &generate_llama_ref);
#endif

#if defined(HAVE_STINGER)
    result.emplace_back("stinger", "Stinger library", &generate_stinger);
    result.emplace_back("stinger-dv", "Stinger with dense vertices", &generate_stinger_dv);
    result.emplace_back("stinger-ref", "Stinger with the GAPBS ref impl.", &generate_stinger_ref);
#endif

#if defined(HAVE_GRAPHONE)
    result.emplace_back("g1.2-cons-sp", "GraphOne, consistency for updates, sparse vertices (vertex dictionary)", &generate_graphone_cons_sp);
    result.emplace_back("g1.2-cons-dv", "GraphOne, consistency for updates, dense vertices", &generate_graphone_cons_dv);
    result.emplace_back("g1.2-bw-sp", "GraphOne, blind writes, sparse vertices (vertex dictionary)", &generate_graphone_bw_sp);
    result.emplace_back("g1.2-bw-dv", "GraphOne, blind writes, dense vertices", &generate_graphone_bw_dv);
    result.emplace_back("g1.2-bw-sp-ignore-build", "GraphOne, blind writes, sparse vertices (vertex dictionary), new deltas/levels cannot be explicitly created", &generate_graphone_bw_sp_ignore_build);
    result.emplace_back("g1.2-bw-dv-ignore-build", "GraphOne, blind writes, dense vertices, new deltas/levels cannot be explicitly created", &generate_graphone_bw_dv_ignore_build);
    result.emplace_back("g1.2-ref", "GraphOne, reference GAP BS for the Graphalytics algorithms", &generate_graphone_ref);
    result.emplace_back("g1.2-ref-ignore-build", "GraphOne, reference GAP BS for the Graphalytics algorithms", &generate_graphone_ref_ignore_build);
#endif

#if defined(HAVE_TESEO)
    result.emplace_back("teseo.2", "Teseo", &generate_teseo);
#endif

    return result;
}

/*****************************************************************************
 *                                                                           *
 *  Base interface                                                           *
 *                                                                           *
 *****************************************************************************/
Interface::Interface(){}
Interface::~Interface(){}
void Interface::on_main_init(int num_threads){ };
void Interface::on_thread_init(int thread_id){ };
void Interface::on_thread_destroy(int thread_id){ } ;
void Interface::on_main_destroy(){ };
bool Interface::has_edge(uint64_t source, uint64_t destination) const {
    return !isnan(get_weight(source, destination));
}
void Interface::dump() const{
    dump_ostream(std::cout);
}
void Interface::dump(const std::string& path) const {
    fstream handle(path.c_str());
    if(!handle.good()) ERROR("[dump] Cannot open the file: `" << path << "'");
    dump_ostream(handle);
    handle.close();
}

void Interface::dump(const char* c_path) const {
    string path = c_path;
    dump(path);
}

bool Interface::is_undirected() const {
    return !is_directed();
}


/*****************************************************************************
 *                                                                           *
 *  Update interface                                                         *
 *                                                                           *
 *****************************************************************************/
bool UpdateInterface::batch(const SingleUpdate* array, size_t array_sz, bool force){
    bool result = true;
    const SingleUpdate* __restrict A = array;
    COUT_DEBUG("batch: " << array_sz << ", force: " << force);

    if(force){

        // now a bit of a hack, we want all updates to succeed. An update may fail if a vertex is still being added
        // by another thread in the meanwhile
        for(uint64_t i = 0; i < array_sz; i++){
            if(A[i].m_weight >= 0){ // insert
                graph::WeightedEdge edge{A[i].m_source, A[i].m_destination, A[i].m_weight};
                if( ! add_edge(edge) ){ // avoid infinite loops/waits
                    auto op = [this](graph::WeightedEdge edge){ return add_edge(edge); };
                    batch_try_again(op, edge);
                }
            } else { // remove
                graph::Edge edge{A[i].m_source, A[i].m_destination};

                if ( ! remove_edge(edge) ){ // avoid infinite loops/waits
                    auto op = [this](graph::Edge edge){ return remove_edge(edge); };
                    batch_try_again(op, edge);
                }
            }
        }
    } else {
        for(uint64_t i = 0; i < array_sz; i++){
            if(A[i].m_weight >= 0){ // insert
                result &= add_edge(graph::WeightedEdge{A[i].m_source, A[i].m_destination, A[i].m_weight});
            } else { // remove
                result &= remove_edge(graph::Edge{A[i].m_source, A[i].m_destination});
            }
        }
    }

    return result;
}

template<typename Action, typename Edge>
void UpdateInterface::batch_try_again(Action action, Edge edge){
    constexpr chrono::seconds timeout = 10min;
    bool result = false;

    auto t_start = chrono::steady_clock::now();
    do {
        result = std::invoke(action, edge);
    } while(!result && chrono::steady_clock::now() - t_start <= timeout);

    if(!result){
        RAISE_EXCEPTION(TimeoutError, "Cannot process the edge " << edge << " after " << timeout.count() << " seconds");
    }
}

void UpdateInterface::load(const string& path) {
    auto reader = reader::Reader::open(path);
    ASSERT(reader->is_directed() == is_directed());
    gfe::graph::WeightedEdge edge;
    while(reader->read(edge)){
        add_vertex(edge.m_source);
        add_vertex(edge.m_destination);
        add_edge(edge);
    }
    build();
}

void UpdateInterface::build(){
    /* nop */
}

uint64_t UpdateInterface::num_levels() const {
    return 0; // by default, we assume that the implementation is not LSM/delta based, and it doesn`t create new levels/deltas/snapshots
}

} // namespace library
