/**
 * The following source code is derived from
 * 1- llama/include/llama/ll_mlcsr_iterator.h
 * 2- llama/benchmark/benchmarks/triangle_counting.h
 * These files report the following attributions:
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

#include <fstream>
#include <iostream>
#include <set>

#include "common/timer.hpp"
#include "utility/timeout_service.hpp"

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
 *  Implementation, directed graph                                           *
 *                                                                           *
 *****************************************************************************/
// this implementation is derived from:
// 1. the class ll_b_triangle_counting_org, as defined in llama/benchmark/benchmarks/triangle_counting.h ; and
// 2. the class ll_common_neighbor_iter, defined in llama/include/llama/ll_mlcsr_iterator.h

namespace { // anonymous
class common_neighbor_iter_directed {
    ll_mlcsr_ro_graph& G;
    const node_t src;
    const node_t dst;

    ll_edge_iterator src_iter;
    bool src_iter_reverse_edge = false;
    ll_edge_iterator dst_iter;
    node_t dst_cur_neighbour = LL_NIL_NODE;
    std::set<node_t> already_visited;

public:

    // graph, source, destination
    common_neighbor_iter_directed(ll_mlcsr_ro_graph& graph, node_t s, node_t d) : G(graph), src(s), dst(d) {
        G.out_iter_begin(src_iter, src);
        G.out_iter_begin(dst_iter, dst);

        edge_t e = G.out_iter_next(dst_iter);
        if(e != LL_NIL_EDGE){
            dst_cur_neighbour = LL_ITER_OUT_NEXT_NODE(G, dst_iter, e);
        }
    }

    node_t get_next() {
        if(dst_cur_neighbour == LL_NIL_NODE){
            // dst is a sink, it has no outgoing edges. No triangles can be built passing through this node
            return LL_NIL_NODE;
        } else if(!src_iter_reverse_edge){
            while(G.out_iter_has_next(src_iter)){
                [[maybe_unused]] edge_t e = G.out_iter_next(src_iter);
                assert(e != LL_NIL_EDGE);
                node_t t = LL_ITER_OUT_NEXT_NODE(G, src_iter, e);
                assert(t != LL_NIL_NODE);
                already_visited.insert(t);
                if (check_common(t)) return t;
            }

            // time to inspect the reverse edges
            G.in_iter_begin_fast(src_iter, src);
            src_iter_reverse_edge = true;
            return get_next();
        } else {
            while(G.in_iter_has_next_fast(src_iter)){
                [[maybe_unused]] edge_t e = G.in_iter_next_fast(src_iter);
                assert(e != LL_NIL_EDGE);
                node_t t = LL_ITER_OUT_NEXT_NODE(G, src_iter, e);
                assert(t != LL_NIL_NODE);
                if(already_visited.count(t)) continue;
                if (check_common(t)) return t;
            }

            return LL_NIL_NODE;
        }
    }

private:
    bool check_common(node_t t) {
        if(t == dst_cur_neighbour) return true;

        node_t candidate = dst_cur_neighbour;

        do {
            while(G.out_iter_has_next(dst_iter)){
                [[maybe_unused]] edge_t e = G.out_iter_next(dst_iter);
                assert(e != LL_NIL_EDGE);
                candidate = LL_ITER_OUT_NEXT_NODE(G, dst_iter, e);
                assert(candidate != LL_NIL_NODE);

                if(candidate == t){ // found
                    dst_cur_neighbour = candidate;

                    return true;
                } else if (candidate == dst_cur_neighbour){ // full iteration => no results
                    return false;
                }
            }

            // restart
            G.out_iter_begin(dst_iter, dst);
        } while (true);

        return false;
    }
};


static
pair<uint64_t, uint64_t> count_triangles_directed(ll_mlcsr_ro_graph& graph, node_t u){
    uint64_t num_triangles = 0;
    set<node_t> neighbours_already_visited;
    uint64_t degree = 0;

    ll_edge_iterator iterator;
    graph.in_iter_begin_fast(iterator, u);
    for (edge_t e = graph.in_iter_next_fast(iterator); e != LL_NIL_EDGE; e = graph.in_iter_next_fast(iterator)) {
        node_t v = LL_ITER_OUT_NEXT_NODE(graph, iterator, e);
        neighbours_already_visited.insert(v);
        degree++;

        common_neighbor_iter_directed w_I(graph, u, v);
        for (node_t w = w_I.get_next(); w != LL_NIL_NODE; w = w_I.get_next()) {
            COUT_DEBUG("Triangle found: " << u << " <- " << v << " -> " << w);
            num_triangles++;
        }
    }

    graph.out_iter_begin(iterator, u);
    for (edge_t e = graph.out_iter_next(iterator); e != LL_NIL_EDGE; e = graph.out_iter_next(iterator)) {
        node_t v = LL_ITER_OUT_NEXT_NODE(graph, iterator, e);
        if(neighbours_already_visited.count(v)) continue; // vertex already processed in the incoming edges iterator
        degree++;

        common_neighbor_iter_directed w_I(graph, u, v);
        for (node_t w = w_I.get_next(); w != LL_NIL_NODE; w = w_I.get_next()) {
            COUT_DEBUG("Triangle found: " << u << " -> " << v << " -> " << w);
            num_triangles++;
        }
    }

    return { num_triangles, degree };
}

static
void llama_execute_lcc_directed(TimeoutService& timer, ll_mlcsr_ro_graph& graph, double* scores){

    #pragma omp for nowait schedule(dynamic,4096)
    for(node_t u = 0; u < graph.max_nodes(); u++){
        if(timer.is_timeout()) continue; // we're done

        if(graph.in_degree(u) + graph.out_degree(u) <= 1){
            scores[u] = 0.0;
        } else {
            auto result = count_triangles_directed(graph, u);
            uint64_t num_triangles = result.first;
            uint64_t degree_no_dups = result.second; // degree of the set of neighbours, ignoring duplicates in the incoming/outgoing edges
            if(degree_no_dups == 1){ // this node has one incoming and one outgoing edge to the same neighbour!
                scores[u] = 0.0;
            } else {
                uint64_t max_num_edges = degree_no_dups * (degree_no_dups -1);
                scores[u] = static_cast<double>(num_triangles) / max_num_edges;
            }
        }
    }
}

} // anonymous namespace

/*****************************************************************************
 *                                                                           *
 *  Implementation, undirected graph                                         *
 *                                                                           *
 *****************************************************************************/
// this implementation is derived from:
// 1. the class ll_b_triangle_counting_org, as defined in llama/benchmark/benchmarks/triangle_counting.h ; and
// 2. the class ll_common_neighbor_iter, defined in llama/include/llama/ll_mlcsr_iterator.h

namespace { // anonymous
class common_neighbor_iter_undirected {
    ll_mlcsr_ro_graph& G;
    const node_t src;
    const node_t dst;

    ll_edge_iterator src_iter;
    bool src_iter_reverse_edge = false;
    ll_edge_iterator dst_iter;
    node_t dst_cur_neighbour = LL_NIL_NODE;
    bool dst_iter_reverse_edge = false;

public:

    // graph, source, destination
    common_neighbor_iter_undirected(ll_mlcsr_ro_graph& graph, node_t s, node_t d) : G(graph), src(s), dst(d) {
        G.out_iter_begin(src_iter, src);
        G.out_iter_begin(dst_iter, dst);

        edge_t e = G.out_iter_next(dst_iter);
        if(e != LL_NIL_EDGE){
            dst_cur_neighbour = LL_ITER_OUT_NEXT_NODE(G, dst_iter, e);
        } else {
            dst_iter_reverse_edge = true;
            G.in_iter_begin_fast(dst_iter, dst);
            G.in_iter_next_fast(dst_iter);

            // assuming that src and dst are neighbours, then dst must have at least one neighbour: src
            e = G.in_iter_has_next_fast(dst_iter);
            assert(e != LL_NIL_EDGE && "dst has no incoming or outgoing edges");
            dst_cur_neighbour = LL_ITER_OUT_NEXT_NODE(G, dst_iter, e);
        }
    }

    node_t get_next() {
        if(!src_iter_reverse_edge){
            while(G.out_iter_has_next(src_iter)){
                [[maybe_unused]] edge_t e = G.out_iter_next(src_iter);
                assert(e != LL_NIL_EDGE);
                node_t t = LL_ITER_OUT_NEXT_NODE(G, src_iter, e);
                assert(t != LL_NIL_NODE);
                if (check_common(t)) return t;
            }

            // time to inspect the reverse edges
            G.in_iter_begin_fast(src_iter, src);
            src_iter_reverse_edge = true;
            return get_next();
        } else {
            while(G.in_iter_has_next_fast(src_iter)){
                [[maybe_unused]] edge_t e = G.in_iter_next_fast(src_iter);
                assert(e != LL_NIL_EDGE);
                node_t t = LL_ITER_OUT_NEXT_NODE(G, src_iter, e);
                assert(t != LL_NIL_NODE);
                if (check_common(t)) return t;
            }

            return LL_NIL_NODE;
        }
    }

private:
    bool check_common(node_t t) {
        if(t == dst_cur_neighbour) return true;

        bool mode_reverse_edge = dst_iter_reverse_edge;
        node_t candidate = dst_cur_neighbour;

        do {
            if(mode_reverse_edge == false){ // outgoing edges
                while(G.out_iter_has_next(dst_iter)){
                    [[maybe_unused]] edge_t e = G.out_iter_next(dst_iter);
                    assert(e != LL_NIL_EDGE);
                    candidate = LL_ITER_OUT_NEXT_NODE(G, dst_iter, e);
                    assert(candidate != LL_NIL_NODE);

                    if(candidate == t){ // found
                        dst_iter_reverse_edge = mode_reverse_edge;
                        dst_cur_neighbour = candidate;

                        return true;
                    } else if (candidate == dst_cur_neighbour){ // full iteration => no results
                        return false;
                    }
                }

                // switch mode
                G.in_iter_begin_fast(dst_iter, dst);
                mode_reverse_edge = true;
            } else {
                while(G.in_iter_has_next_fast(dst_iter)){
                    [[maybe_unused]] edge_t e = G.in_iter_next_fast(dst_iter);
                    assert(e != LL_NIL_EDGE);
                    candidate = LL_ITER_OUT_NEXT_NODE(G, dst_iter, e);
                    assert(candidate != LL_NIL_NODE);

                    if(candidate == t){ // found
                        dst_iter_reverse_edge = mode_reverse_edge;
                        dst_cur_neighbour = candidate;

                        return true;
                    } else if (candidate == dst_cur_neighbour){ // full iteration => no results
                        return false;
                    }
                }

                // switch mode
                G.out_iter_begin(dst_iter, dst);
                mode_reverse_edge = false;
            }
        } while (true);

        return false;
    }
};


static
uint64_t count_triangles_undirected(ll_mlcsr_ro_graph& graph, node_t u){
    uint64_t result = 0;

    ll_edge_iterator iterator;
    graph.in_iter_begin_fast(iterator, u);
    for (edge_t e = graph.in_iter_next_fast(iterator); e != LL_NIL_EDGE; e = graph.in_iter_next_fast(iterator)) {
        node_t v = LL_ITER_OUT_NEXT_NODE(graph, iterator, e);
        common_neighbor_iter_undirected w_I(graph, u, v);
        for (node_t w = w_I.get_next(); w != LL_NIL_NODE; w = w_I.get_next()) {
            COUT_DEBUG("Triangle found: " << u << " - " << v << " - " << w);
            result++;
        }
    }

    graph.out_iter_begin(iterator, u);
    for (edge_t e = graph.out_iter_next(iterator); e != LL_NIL_EDGE; e = graph.out_iter_next(iterator)) {
        node_t v = LL_ITER_OUT_NEXT_NODE(graph, iterator, e);
        common_neighbor_iter_undirected w_I(graph, u, v);
        for (node_t w = w_I.get_next(); w != LL_NIL_NODE; w = w_I.get_next()) {
            COUT_DEBUG("Triangle found: " << u << " - " << v << " - " << w);
            result++;
        }
    }

    return result;
}

static
void llama_execute_lcc_undirected(TimeoutService& timer, ll_mlcsr_ro_graph& graph, double* scores){

    #pragma omp for nowait schedule(dynamic,4096)
    for(node_t u = 0; u < graph.max_nodes(); u++){
        if(timer.is_timeout()) continue; // we're done

        uint64_t degree = graph.in_degree(u) + graph.out_degree(u);
        if(degree <= 1){
            scores[u] = 0.0;
        } else {
            uint64_t num_triangles = count_triangles_undirected(graph, u);
            uint64_t max_num_edges = degree * (degree -1);
            scores[u] = static_cast<double>(num_triangles) / max_num_edges;
        }
    }

}

} // anonymous namespace

/*****************************************************************************
 *                                                                           *
 *  Wrapper                                                                  *
 *                                                                           *
 *****************************************************************************/
namespace library {

void LLAMAClass::lcc(const char* dump2file){
    TimeoutService timeout { m_timeout };
    Timer timer; timer.start();

    // retrieve the latest snapshot the internal source_vertex_id
    shared_lock<shared_mutex> slock(m_lock_checkpoint);
    auto graph = get_snapshot();
//    dump_snapshot(graph);
    slock.unlock();

    // execute the LCC algorithm
    unique_ptr<double[]> ptr_scores { new double[graph.max_nodes()] };
    double* scores = ptr_scores.get();
    if(is_directed()){
        llama_execute_lcc_directed(timeout, graph, /* output */ scores);
    } else {
        llama_execute_lcc_undirected(timeout, graph, /* output */ scores);
    }

    if(timeout.is_timeout()){
        RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);
    }

    // translate from llama vertex ids to external vertex ids
    auto names = graph.get_node_property_64(g_llama_property_names);
    assert(names != nullptr && "Wrong string ID to refer the property attached to the vertices");
    cuckoohash_map</* external id */ uint64_t, /* score */ double> external_ids;
    #pragma omp parallel for
    for(node_t llama_node_id = 0; llama_node_id < graph.max_nodes(); llama_node_id++){
        // first, does this node exist (or it's a gap?)
        // this is a bit of a stretch: the impl~ from llama assumes that a node does not exist only if it does not have any incoming or outgoing edges.
        if(!graph.node_exists(llama_node_id)) continue;

        // second, what's it's real node ID, in the external domain (e.g. user id)
        uint64_t external_node_id = names->get(llama_node_id);

        // third, its score
        double score = scores[llama_node_id];

        // finally, register the association
        external_ids.insert(external_node_id, score);
    }

    if(timeout.is_timeout()){
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

