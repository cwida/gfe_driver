/**
 * The following source code is derived from
 * 1- llama/include/llama/ll_dfs_template.h
 * 2- llama/benchmark/benchmarks/tarjan_scc.h
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

/*
 * This file was adapted from Green-Marl, which includes the following notice:
 *
 * Copyright (c) 2011-2012 Stanford University, unless otherwise specified.
 * All rights reserved.
 *
 * This software was developed by the Pervasive Parallelism Laboratory of
 * Stanford University, California, USA.
 *
 * Permission to use, copy, modify, and distribute this software in source or
 * binary form for any purpose with or without fee is hereby granted, provided
 * that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of Stanford University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
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

#include <cstring>
#include <fstream>
#include <iostream>
#include <set>
#include <shared_mutex> // shared_lock
#include <unordered_set>
#include <utility>
#include <vector>

#include "llama/ll_seq.h"

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
 *  LLAMA's implementation                                                   *
 *                                                                           *
 *****************************************************************************/
namespace { // anonymous

template<class Graph, bool has_pre_visit, bool has_post_visit, bool has_navigator>
class ll_dfs_template {
protected:
    gfe::utility::TimeoutService& timeout_srv;
    node_t root;
    Graph& G;

    // stack implementation
    node_t stack_ptr;
    std::vector<std::pair<ll_edge_iterator, bool>> stack;
    ll_edge_iterator iter;

    // visited set implementation
    node_t cnt;
    unsigned char* visited_bitmap;
    std::unordered_set<node_t> visited_small;
    bool is_small;
    int THRESHOLD_LARGE;
    static const node_t INVALID_NODE = -1;

    bool use_reverse_edge = false; // whether we are visiting the incoming edges in the current node being iterated

protected:
    virtual void visit_pre(node_t t) { };
    virtual void visit_post(node_t t) { };
    virtual bool check_navigator(node_t t, edge_t idx) { return true; };

public:
    ll_dfs_template(gfe::utility::TimeoutService& timeout_srv, Graph& _G) : timeout_srv(timeout_srv), G(_G) {
        visited_bitmap = NULL; // bitmap
    }

    virtual ~ll_dfs_template() {
        delete visited_bitmap;
    }

    void prepare(node_t root_node) {
        root = root_node;
        cnt = 0;
        visited_small.clear();

        is_small = true;;
        use_reverse_edge = true;
        iter.node = INVALID_NODE;
        THRESHOLD_LARGE = std::max((int)(G.max_nodes()*0.1), 4096);
    }

    void do_dfs() {
        COUT_DEBUG("init: " << root);
        enter_node(root);
        main_loop();
    }

private:
    void enter_node(node_t n) {
        COUT_DEBUG("[" << stack.size() << "] node: " << n);

        // push current node
        stack.push_back({ iter, use_reverse_edge });

        G.out_iter_begin(iter, n);
        use_reverse_edge = false;

        // mark visited
        add_visited(n);
        cnt++;
        if (cnt == THRESHOLD_LARGE) {
            prepare_large();
        }

        if (has_pre_visit) visit_pre(n);
    }


    void prepare_large() {
        delete[] visited_bitmap;

        visited_bitmap = new unsigned char[(G.max_nodes() + 7) / 8];

        #pragma omp parallel for schedule(dynamic,16384)
        for (int i = 0; i < (G.max_nodes() + 7) / 8; i++)
            visited_bitmap[i] = 0;

        std::unordered_set<node_t>::iterator I;
        for (I = visited_small.begin(); I != visited_small.end(); I++) {
            node_t u = *I;
            _ll_set_bit(visited_bitmap, u);
        }
        is_small = false;
        stack.reserve(G.max_nodes());
    }


    void exit_node(node_t n) {
        if (has_post_visit) visit_post(n);
        COUT_DEBUG("[" << stack.size() -1 << "] node: " << n);
        auto& pair = stack.back();
        iter = pair.first;
        use_reverse_edge = pair.second;
        stack.pop_back();
    }

    void main_loop() {
        //----------------------------------
        // Repeat until stack is empty
        //----------------------------------
        while (iter.node != INVALID_NODE && !timeout_srv.is_timeout()) {
            //----------------------------------
            // Every neighbor has been visited
            //----------------------------------
            if (iter.edge == LL_NIL_EDGE) {
                if(use_reverse_edge == false){
                    COUT_DEBUG("[" << stack.size() -1 << "] switch: " << iter.node);
                    G.in_iter_begin_fast(iter, iter.node);
                    use_reverse_edge = true;
                } else { // we are done visiting iter.node
                    exit_node(iter.node);
                }

            } else {
                //----------------------------------
                // check every non-visited neighbor
                //----------------------------------
                node_t z;
                edge_t e;
                if (use_reverse_edge) {
                    e = G.in_iter_next_fast(iter);
                    COUT_DEBUG("[" << stack.size() -1 << "] in: " << iter.node << " <- " << iter.last_node);
                } else {
                    e = G.out_iter_next(iter);
                    COUT_DEBUG("[" << stack.size() -1 << "] out: " << iter.node << " -> " << iter.last_node);
                }
                assert(e != LL_NIL_EDGE);
                z = iter.last_node;

                if (has_visited(z)) {
                    continue;
                }

                if (has_navigator) {
                    if (check_navigator(z, e) == false) {
                        continue;
                    }
                }
                enter_node(z);
                continue;
            }
        }
    }

    void add_visited(node_t n) {
        if (is_small)
            visited_small.insert(n);
        else
            _ll_set_bit(visited_bitmap, n);
    }

    bool has_visited(node_t n) {
        if (is_small) {
            return (visited_small.find(n) != visited_small.end());
        } else {
            return _ll_get_bit(visited_bitmap, n);
        }
    }
};


template <class Graph>
class WeaklyConnectedComponents : public ll_dfs_template <Graph, true, false, false> {
public:
    WeaklyConnectedComponents(TimeoutService& timeout_srv, Graph& _G, node_t*& _G_SCC, node_t _n)
        : ll_dfs_template<Graph, true, false, false>(timeout_srv, _G),  G(_G), G_SCC(_G_SCC), n(_n){

    }

private:  // list of variables
    Graph& G;
    node_t*& G_SCC;
    node_t n;

protected:
    virtual void visit_pre(node_t t) {
        G_SCC[t] = n;
    }
};

static
void llama_execute_wcc(TimeoutService& timeout_srv, ll_mlcsr_ro_graph& graph, /* output array, already allocated */ node_t* G_SCC){

    #pragma omp parallel for
    for (node_t t0 = 0; t0 < graph.max_nodes(); t0 ++) {
        graph.set_node_prop(G_SCC, t0, LL_NIL_NODE);
    }

    for (node_t n = 0; n < graph.max_nodes(); n++ ){
        if (G_SCC[n] == LL_NIL_NODE) {
            WeaklyConnectedComponents<ll_mlcsr_ro_graph> _DFS(timeout_srv, graph, G_SCC,  n);
            _DFS.prepare(n);
            _DFS.do_dfs();
        }
    }
}

} // anonymous namespace

/*****************************************************************************
 *                                                                           *
 *  Wrapper                                                                  *
 *                                                                           *
 *****************************************************************************/
namespace gfe::library {

void LLAMAClass::wcc(const char* dump2file){
    utility::TimeoutService timeout_srv { m_timeout };
    Timer timer; timer.start();

    // retrieve the latest snapshot the internal source_vertex_id
    shared_lock<shared_mutex_t> slock(m_lock_checkpoint);
    auto graph = get_snapshot();
//    dump_snapshot(graph);
    slock.unlock();

    // execute the WCC algorithm
    unique_ptr<node_t[]> ptr_components { new node_t[graph.max_nodes()] };
    node_t* components = ptr_components.get();
    llama_execute_wcc(timeout_srv, graph, /* output */ components);
    if(timeout_srv.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    // translate from llama ids to external vertex ids
    auto external_ids = translate(graph, components);
    if(timeout_srv.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    if(dump2file != nullptr) // store the results in the given file
        save_results(external_ids, dump2file);
}

} // namespace

