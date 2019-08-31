/**
 * The following source code is derived from
 * 1- llama/benchmark/benchmarks/sssp.h
 * This file reports the following attributions:
 */

/*
 * LLAMA Graph Analytics
 *
 * Copyright 2014
 *      The President and Fellows of Harvard College.
 *
 * Copyright 2014
 *      Oracle Labs.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE UNIVERSITY OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

#include "llama_class.hpp"
#include "llama_internal.hpp"

#include "common/timer.hpp"
#include "utility/timeout_service.hpp"

#include <fstream>
#include <memory>
#include <iostream>

using namespace common;
using namespace std;
using namespace utility;

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
 *  SSSP                                                                     *
 *                                                                           *
 *****************************************************************************/
// the following source code is based on the class `ll_b_sssp_weighted', located
// in llama/benchmark/benchmarks/sssp.h
// @param G_dist is the output array with the distances, it must be already allocated (but not initialised) and its length equal at least to graph.max_nodes()
static void llama_execute_sssp(TimeoutService& timer, ll_mlcsr_ro_graph& graph, node_t root, const char* weights_property_name, bool is_graph_undirected, double* __restrict G_dist){
    ll_mlcsr_edge_property<double>& G_len = *reinterpret_cast<ll_mlcsr_edge_property<double>*>(graph.get_edge_property_64(weights_property_name));
    ll_spinlock_table lt;
    ll_memory_helper m;

    bool fin = false ;

    bool* G_updated = m.allocate<bool>(graph.max_nodes());
    bool* G_updated_nxt = m.allocate<bool>(graph.max_nodes());
    unique_ptr<double[]> ptr_G_dist_nxt { new double[graph.max_nodes()] }; // avoid memory leaks
    double* G_dist_nxt = ptr_G_dist_nxt.get();

    fin = false;

    // Initialisation
    #pragma omp parallel for
    for (node_t t0 = 0; t0 < graph.max_nodes(); t0++) {
        G_dist[t0] = (t0 == root) ? 0.0 : numeric_limits<double>::infinity();
        G_updated[t0] = (t0 == root) ? true : false;
        G_dist_nxt[t0] = G_dist[t0];
        G_updated_nxt[t0] = G_updated[t0];
    }

    // Graph visit
    while (!fin && !timer.is_timeout()) {
        COUT_DEBUG("init iteration");
        bool __E8 = false; // no idea what this is

        fin = true;
        __E8 = false;

        #pragma omp parallel for schedule(dynamic,4096)
        for (node_t n = 0; n < graph.max_nodes(); n++) {
            if (G_updated[n]) {
                ll_edge_iterator iter;
                graph.out_iter_begin(iter, n);
                for (edge_t s_idx = graph.out_iter_next(iter); s_idx != LL_NIL_EDGE; s_idx = graph.out_iter_next(iter)) {
                    node_t s = LL_ITER_OUT_NEXT_NODE(graph, iter, s_idx);
                    edge_t e;

                    e = s_idx;
                    { // argmin(argmax) - test and test-and-set
                        double G_dist_nxt_new = G_dist[n] + G_len[e];
                        if (G_dist_nxt[s] > G_dist_nxt_new) {
                            bool G_updated_nxt_arg = true;
                            lt.acquire_for(s);
                            if (G_dist_nxt[s] > G_dist_nxt_new) {
                                COUT_DEBUG("out edge: " << n << " -> (" << G_dist_nxt_new << ") -> " << s);
                                G_dist_nxt[s] = G_dist_nxt_new;
                                G_updated_nxt[s] = G_updated_nxt_arg;
                            }
                            lt.release_for(s);
                        }
                    }
                }

                // in undirected graphs, the assumption is that both llama's incoming & outgoing edges form the set of undirected edges of a node
                if(is_graph_undirected){ // as above
                    graph.in_iter_begin_fast(iter, n);
                    for (edge_t s_idx = graph.in_iter_next_fast(iter); s_idx != LL_NIL_EDGE;  s_idx = graph.in_iter_next_fast(iter)) {
                        node_t s = LL_ITER_OUT_NEXT_NODE(graph, iter, s_idx);
                        edge_t e;

                        e = graph.in_to_out(s_idx); // the weight is only a property of the outgoing edges
                        { // argmin(argmax) - test and test-and-set

                            double G_dist_nxt_new = G_dist[n] + G_len[e];
                            if (G_dist_nxt[s] > G_dist_nxt_new) {
                                bool G_updated_nxt_arg = true;
                                lt.acquire_for(s);
                                if (G_dist_nxt[s] > G_dist_nxt_new) {
                                    COUT_DEBUG("in edge: " << n << " -> (" << G_dist_nxt_new << ") -> " << s);
                                    G_dist_nxt[s] = G_dist_nxt_new;
                                    G_updated_nxt[s] = G_updated_nxt_arg;
                                }
                                lt.release_for(s);
                            }
                        }
                    }
                } // end if (is_graph_undirected)
            }
        }


        #pragma omp parallel
        {
            bool __E8_prv = false ;

            __E8_prv = false ;

            #pragma omp for nowait
            for (node_t t4 = 0; t4 < graph.max_nodes(); t4 ++){
                G_dist[t4] = G_dist_nxt[t4];
                G_updated[t4] = G_updated_nxt[t4];
                G_updated_nxt[t4] = false;
                __E8_prv = __E8_prv || G_updated[t4] ;
            }
            ATOMIC_OR(&__E8, __E8_prv);
        }
        fin =  !__E8 ;
    }
}


/*****************************************************************************
 *                                                                           *
 *  SSSP                                                                     *
 *                                                                           *
 *****************************************************************************/

namespace library {
void LLAMAClass::sssp(uint64_t external_source_vertex_id, const char* dump2file){
    utility::TimeoutService timeout_srv { m_timeout };
    Timer timer; timer.start();

    // retrieve the latest snapshot the internal source_vertex_id
    shared_lock<shared_mutex> slock(m_lock_checkpoint);
    auto graph = get_snapshot();
//    dump_snapshot(graph);
    int64_t llama_source_vertex_id;
    if (! m_vmap_read_only.find(external_source_vertex_id, llama_source_vertex_id) ){ // side effect, it assigns llama_source_vertex_id
        ERROR("The given vertex does not exist: " << external_source_vertex_id);
    }
    slock.unlock();

    // execute the SSSP algorithm
    unique_ptr<double[]> ptr_distances { new double[graph.max_nodes()] };
    double* distances = ptr_distances.get();
    llama_execute_sssp(timeout_srv, graph, llama_source_vertex_id, g_llama_property_weights, is_undirected(), /* output */ distances);

    if(timeout_srv.is_timeout()){
        RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);
    }

    // translate from llama ids to external vertex ids
    auto names = graph.get_node_property_64(g_llama_property_names);
    assert(names != nullptr && "Wrong string ID to refer the property attached to the vertices");
    cuckoohash_map</* external id */ uint64_t, /* distance */ double> external_ids;
    #pragma omp parallel for
    for(node_t llama_node_id = 0; llama_node_id < graph.max_nodes(); llama_node_id++){
        // first, does this node exist (or it's a gap?)
        // this is a bit of a stretch: the impl~ from llama assumes that a node does not exist only if it does not have any incoming or outgoing edges.
        if(!graph.node_exists(llama_node_id)) continue;

        // second, what's it's real node ID, in the external domain (e.g. user id)
        uint64_t external_node_id = names->get(llama_node_id);

        // third, its distance
        double distance = distances[llama_node_id];

        // finally, register the association
        external_ids.insert(external_node_id, distance);
    }

    if(timeout_srv.is_timeout()){
        RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);
    }

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

} // namespace

