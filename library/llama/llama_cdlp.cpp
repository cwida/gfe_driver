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

#include <fstream>
#include <iostream>
#include <shared_mutex> // shared_lock
#include <unordered_map>

#include "common/timer.hpp"
#include "third-party/libcuckoo/cuckoohash_map.hh"
#include "utility/timeout_service.hpp"

using namespace common;
using namespace gfe::utility;
using namespace libcuckoo;
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
 *  CDLP                                                                     *
 *                                                                           *
 *****************************************************************************/
namespace gfe::library {

void LLAMAClass::cdlp(uint64_t max_iterations, const char* dump2file){
#if defined(LL_COUNTERS)
    ll_clear_counters();
#endif

    TimeoutService timeout { m_timeout };
    Timer timer; timer.start();

    // retrieve the latest snapshot the internal source_vertex_id
    shared_lock<shared_mutex_t> slock(m_lock_checkpoint);
    auto graph = get_snapshot();
//    dump_snapshot(graph);
    slock.unlock(); // here we lose the ability to refer to m_vmap_read_only from now on

    // execute the CDLP algortihm
    unique_ptr<uint64_t[]> ptr_labels = cdlp_impl(timeout, graph, max_iterations);
    uint64_t* labels = ptr_labels.get();

    if(timeout.is_timeout()){
        RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);
    }

    // translate from llama vertex ids to external vertex ids
    auto names = graph.get_node_property_64(g_llama_property_names);
    assert(names != nullptr && "Wrong string ID to refer the property attached to the vertices");
    cuckoohash_map</* external id */ uint64_t, /* score */ uint64_t> external_ids;
    #pragma omp parallel for
    for(node_t llama_node_id = 0; llama_node_id < graph.max_nodes(); llama_node_id++){
        // first, does this node exist (or it's a gap?)
        // this is a bit of a stretch: the impl~ from llama assumes that a node does not exist only if it does not have any incoming or outgoing edges.
        if(!graph.node_exists(llama_node_id)) continue;

        // second, what's it's real node ID, in the external domain (e.g. user id)
        uint64_t external_node_id = names->get(llama_node_id);

        // third, its label
        uint64_t label = labels[llama_node_id];

        // finally, register the association
        external_ids.insert(external_node_id, label);
    }

    if(timeout.is_timeout()){
        RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);
    }

#if defined(LL_COUNTERS)
    ll_print_counters(stdout);
#endif

    // store the results in the given file
    if(dump2file != nullptr){
        COUT_DEBUG("save the results to: " << dump2file)
        fstream handle(dump2file, ios_base::out);
        if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

        auto hashtable = external_ids.lock_table();

        for(const auto& keyvaluepair : hashtable){
            handle << keyvaluepair.first << " " << keyvaluepair.second << "\n";
        }

        handle.close();
    }
}

unique_ptr<uint64_t[]> LLAMAClass::cdlp_impl(TimeoutService& timer, ll_mlcsr_ro_graph& graph, uint64_t max_iterations){
    unique_ptr<uint64_t[]> ptr_labels0 { new uint64_t[graph.max_nodes()] };
    unique_ptr<uint64_t[]> ptr_labels1 { new uint64_t[graph.max_nodes()] };
    uint64_t* labels0 = ptr_labels0.get(); // current labels
    uint64_t* labels1 = ptr_labels1.get(); // labels for the next iteration

    // initialisation
    auto names = graph.get_node_property_64(g_llama_property_names);
    #pragma omp parallel for
    for(node_t n = 0; n < graph.max_nodes(); n++){
        labels0[n] = names->get(n);
    }

    // algorithm pass
    bool change = true;
    uint64_t current_iteration = 0;
    while(current_iteration < max_iterations && change && !timer.is_timeout()){
        change = false; // reset the flag

        #pragma omp parallel for shared(change) schedule(dynamic, 64)
        for(node_t node_id = 0; node_id < graph.max_nodes(); node_id++){

            unordered_map<int64_t, uint64_t> histogram;
            ll_edge_iterator iterator;

            // compute the histogram from both the outgoing & incoming edges. The aim is to find the number of each label is shared among
            // the neighbours of node_id
            graph.out_iter_begin(iterator, node_id);
            for (edge_t e = graph.out_iter_next(iterator); e != LL_NIL_EDGE; e = graph.out_iter_next(iterator)) {
                node_t neighbour_id = LL_ITER_OUT_NEXT_NODE(graph, iterator, e);
                histogram[labels0[neighbour_id]]++;

            }

            // cfr. Spec v0.9 pp 14 "If the graph is directed and a neighbor is reachable via both an incoming and
            // outgoing edge, its label will be counted twice"
            graph.in_iter_begin_fast(iterator, node_id);
            for (edge_t e = graph.in_iter_next_fast(iterator); e != LL_NIL_EDGE; e = graph.in_iter_next_fast(iterator)) {
                node_t neighbour_id = LL_ITER_OUT_NEXT_NODE(graph, iterator, e);
                histogram[labels0[neighbour_id]]++;
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

            labels1[node_id] = label_max;
            change |= (labels0[node_id] != labels1[node_id]);
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

} // namespace

