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
namespace gfe{ extern mutex _log_mutex [[maybe_unused]]; }
#define COUT_DEBUG_FORCE(msg) { std::scoped_lock<std::mutex> lock{::gfe::_log_mutex}; std::cout << "[LLAMA::" << __FUNCTION__ << "] " << msg << std::endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif

namespace gfe::library {

/*****************************************************************************
 *                                                                           *
 *  Pagerank                                                                 *
 *                                                                           *
 *****************************************************************************/

void LLAMAClass::pagerank(uint64_t num_iterations, double damping_factor, const char* dump2file){
    TimeoutService timeout_srv { m_timeout };
    Timer timer; timer.start();

    // retrieve the latest snapshot
    shared_lock<shared_mutex_t> slock(m_lock_checkpoint);
    auto graph = get_snapshot();
//    dump_snapshot(graph);
    auto current_num_vertices = num_vertices();
    slock.unlock();

    // execute the Pagerank algorithm
    unique_ptr<double[]> ptr_rank { new double[graph.max_nodes()] };
    double* rank = ptr_rank.get();
    pagerank_impl(timeout_srv, graph, current_num_vertices, num_iterations, damping_factor, /* output */ rank);
    if(timeout_srv.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // translate from llama ids to external vertex ids
    auto external_ids = translate(graph, rank);
    if(timeout_srv.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    if(dump2file != nullptr) // store the results in the given file
        save_results(external_ids, dump2file);
}

// Implementation derived from llama/benchmark/benchmarks/pagerank.h, class ll_b_pagerank_pull_ext
void LLAMAClass::pagerank_impl(utility::TimeoutService& timer, ll_mlcsr_ro_graph& graph, uint64_t num_vertices, uint64_t num_iterations, double d, /* output array, already allocated */ double* G_pg_rank){
    COUT_DEBUG("num_vertices: " << num_vertices << ", num_iterations: " << num_iterations << ", damping sum: " << d);

    ll_memory_helper m; // memory pool, it releases the allocated resources at dtor

    double* G_pg_rank_nxt = m.allocate<double>(graph.max_nodes());
    double N = num_vertices;

    // init
    #pragma omp parallel for
    for (node_t t0 = 0; t0 < graph.max_nodes(); t0++) {
        G_pg_rank[t0] = 1.0 / N;
        G_pg_rank_nxt[t0] = 0.0;
    }

    // iterations
    uint64_t current_iteration = 0;
    while(current_iteration < num_iterations && !timer.is_timeout()){

        double dsum = 0.0; // dangling sum
        #pragma omp parallel for reduction(+:dsum)
        for(node_t t = 0; t < graph.max_nodes(); t ++){
            // okay, LLAMA does not know if a vertex does not exist or it's just a gap / empty cell in its arrays
            // we assume all nodes with no incoming or outgoing edges are gaps
            if(graph.out_degree(t) + graph.in_degree(t) == 0) continue;

            if(get_read_store_outdegree(graph, t) == 0){ // this is a sink
                dsum += G_pg_rank[t];
            }

            double score = 0;
            ll_edge_iterator iterator;
            graph.in_iter_begin_fast(iterator, t);
            for (edge_t e = graph.in_iter_next_fast(iterator); e != LL_NIL_EDGE; e = graph.in_iter_next_fast(iterator)) {
                node_t n = LL_ITER_OUT_NEXT_NODE(*G, iterator, e);
                score += G_pg_rank[n] / get_read_store_outdegree(graph, n);
            }

            if(!m_is_directed){ // the actual set of edges in undirected graphs is composed by both LLAMA's incoming and outgoing edges
                graph.out_iter_begin(iterator, t);
                for (edge_t e = graph.out_iter_next(iterator); e != LL_NIL_EDGE; e = graph.out_iter_next(iterator)) {
                    node_t n = LL_ITER_OUT_NEXT_NODE(*G, iterator, e);
                    score += G_pg_rank[n] / get_read_store_outdegree(graph, n);
                }
            }

            G_pg_rank_nxt[t] = score;
        }

        // rank <- rank_next; rank_next <- 0;
        #pragma omp for schedule(dynamic,4096)
        for (node_t t = 0; t < graph.max_nodes(); t ++) {
            G_pg_rank[t] = (1-d) / N + d * G_pg_rank_nxt[t] + /* dangling sum */ d * dsum / N;
            G_pg_rank_nxt[t] = 0.0;
        }

        // next iteration
        current_iteration++;
    }
}

} // namespace

