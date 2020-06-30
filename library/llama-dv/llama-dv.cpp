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
 *
 * Part of the following source code from the LLAMA's source code, which includes
 * the following copyright notices:
 *
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
 *
 * ---------------------------------------------------------
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

#include "llama-dv.hpp"
#include "../llama/llama_internal.hpp"

#include <cstdlib> // exit, debug only
#include <fstream>
#include <iostream>
#include <mutex>
#include <set>
#include <shared_mutex> // shared_lock

#include "common/timer.hpp"
#include "utility/timeout_service.hpp"

using namespace common;
using namespace std;

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG
namespace gfe { extern mutex _log_mutex [[maybe_unused]]; }
#define COUT_DEBUG_FORCE(msg) { std::scoped_lock<std::mutex> lock{gfe::_log_mutex}; std::cout << "[LLAMA_DV::" << __FUNCTION__ << "] " << msg << std::endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif

// External counters, to profile ll_writable_graph#add_edge_if_not_exists and #build
// defined in llama_internal.hpp
#if defined(LL_PROFILE_UPDATES)
extern common::Timer<true> g_llama_total_time;
extern std::atomic<uint64_t> g_llama_add_edge_total_nanosecs;
extern std::atomic<uint64_t> g_llama_build_nanosecs;
extern void llama_add_edge_print_stats();
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

namespace gfe::library {
/*****************************************************************************
 *                                                                           *
 *  Initialisation                                                           *
 *                                                                           *
 *****************************************************************************/

LLAMA_DV::LLAMA_DV(bool is_directed, bool blind_writes) : m_is_directed(is_directed), m_blind_writes(blind_writes) {
    m_db = new ll_database();

    auto& csr = m_db->graph()->ro_graph();
    csr.create_uninitialized_edge_property_64(g_llama_property_weights, LL_T_DOUBLE);
}


LLAMA_DV::~LLAMA_DV(){
#if defined(LL_PROFILE_UPDATES)
    llama_add_edge_print_stats();
#endif

    delete m_db; m_db = nullptr;
}

/******************************************************************************
 *                                                                            *
 *  Properties                                                                *
 *                                                                            *
 *****************************************************************************/
uint64_t LLAMA_DV::num_edges() const {
    shared_lock<shared_mutex_t> cplock(m_lock_checkpoint); // forbid any checkpoint now
    return m_num_edges;
}

uint64_t LLAMA_DV::num_vertices() const {
    shared_lock<shared_mutex_t> cplock(m_lock_checkpoint); // forbid any checkpoint now
    return m_db->graph()->max_nodes();
}

uint64_t LLAMA_DV::num_levels() const{
    shared_lock<shared_mutex_t> cplock(m_lock_checkpoint); // forbid any checkpoint now
    return m_db->graph()->ro_graph().num_levels();
}

bool LLAMA_DV::is_directed() const {
    return m_is_directed;
}

bool LLAMA_DV::has_vertex(uint64_t vertex_id) const {
    return m_db->graph()->node_exists(vertex_id);
}

double LLAMA_DV::get_weight(uint64_t source, uint64_t destination) const {
    if(!m_is_directed && source > destination){ swap(source, destination); }

    edge_t edge_id = m_db->graph()->find(source, source);
    if(edge_id == LL_NIL_EDGE) return numeric_limits<double>::quiet_NaN(); // the edge does not exist
    return get_out_edge_weight(* (m_db->graph()), edge_id);
}

void LLAMA_DV::set_timeout(uint64_t seconds) {
    m_timeout = chrono::seconds{ seconds };
}

uint64_t LLAMA_DV::get_read_store_outdegree(ll_mlcsr_ro_graph& snapshot, int64_t vertex_id) const{
    if(m_is_directed){
        return snapshot.out_degree(vertex_id);
    } else {
        return snapshot.out_degree(vertex_id) + snapshot.in_degree(vertex_id);
    }
}

uint64_t LLAMA_DV::get_write_store_outdegree(int64_t vertex_id) const{
    ll_writable_graph* graph = m_db->graph();

    if(m_is_directed){
        return graph->out_degree(vertex_id);
    } else {
        return graph->out_degree(vertex_id) + graph->in_degree(vertex_id);
    }
}

ll_mlcsr_ro_graph LLAMA_DV::get_snapshot() const {
    if(num_levels() == 0)
        ERROR("There are no levels/deltas/snapshots available. Create them with #build()");

    return ll_mlcsr_ro_graph{ &(m_db->graph()->ro_graph()) , static_cast<int>(num_levels()) -1 };
}

/******************************************************************************
 *                                                                            *
 *  Updates                                                                   *
 *                                                                            *
 *****************************************************************************/
bool LLAMA_DV::add_vertex(uint64_t vertex_id){
    shared_lock<shared_mutex_t> cplock(m_lock_checkpoint); // forbid any checkpoint now
    COUT_DEBUG("vertex_id: " << vertex_id);

    return m_db->graph()->add_node(vertex_id);
}

bool LLAMA_DV::remove_vertex(uint64_t vertex_id){
    shared_lock<shared_mutex_t> cplock(m_lock_checkpoint); // forbid any checkpoint now
    COUT_DEBUG("vertex_id: " << vertex_id);

    m_db->graph()->delete_node(vertex_id);
    return true; // lie
}

bool LLAMA_DV::add_edge(graph::WeightedEdge e){
    shared_lock<shared_mutex_t> cplock(m_lock_checkpoint); // forbid any checkpoint now
    COUT_DEBUG("edge: " << e);
    return add_edge0(e.source(), e.destination(), e.weight());
}

bool LLAMA_DV::add_edge_v2(gfe::graph::WeightedEdge e){
    shared_lock<shared_mutex_t> cplock(m_lock_checkpoint); // forbid any checkpoint now
    COUT_DEBUG("edge: " << e);

    node_t source = (node_t) e.source();
    m_db->graph()->add_node(source);
    node_t destination = (node_t) e.destination();
    m_db->graph()->add_node(destination);

    return add_edge0(source, destination, e.weight());
}

bool LLAMA_DV::add_edge0(int64_t source, int64_t destination, double weight){
#if defined(LL_PROFILE_UPDATES)
    auto t_start = chrono::steady_clock::now();
#endif

    if(!m_is_directed && source > destination){
        std::swap(source, destination);
    }

    edge_t edge_id;
    bool inserted = true;

    if(m_blind_writes){ // blind write, assume that an edge source -> destination does not already exist
        edge_id = m_db->graph()->add_edge(source, destination);
    } else {
        inserted = m_db->graph()->add_edge_if_not_exists(source, destination, &edge_id);
    }

    // thread unsafe, this should really still be under the same latch of add_edge_if_not_exists
    if(inserted){
        m_db->graph()->get_edge_property_64(g_llama_property_weights)->set(edge_id, *reinterpret_cast<uint64_t*>(&(weight)));
    }

#if defined(LL_PROFILE_UPDATES)
    common::compiler_barrier();
    g_llama_add_edge_total_nanosecs += chrono::duration_cast<chrono::nanoseconds>(chrono::steady_clock::now() - t_start).count();
#endif

    return inserted;
}


bool LLAMA_DV::remove_edge(graph::Edge e){
    COUT_DEBUG("edge: " << e);
    shared_lock<shared_mutex_t> cplock(m_lock_checkpoint); // forbid any checkpoint now

    uint64_t source = e.source();
    uint64_t destination = e.destination();
    if(!m_is_directed && source > destination){
        std::swap(source, destination);
    }

    /*
     * m_db->graph()->delete_edge(source, edge_id) requires the id of the edge to remove, which we don't know. The
     * sequence m_db->graph()->delete_edge( source, m_db->graph()->find(source, destination) ) is not thread safe, it
     * would still demand a lock to be correct. Therefore use m_db->graph()->delete_edge_if_exists(...)
     */
    return m_db->graph()->delete_edge_if_exists(source, destination);
}

void LLAMA_DV::build(){
    scoped_lock<shared_mutex_t> xlock(m_lock_checkpoint);
#if defined(LL_PROFILE_UPDATES)
    auto t_start = chrono::steady_clock::now();
#endif
    COUT_DEBUG("build");

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
 *  Overhead to create new delta levels                                       *
 *                                                                           *
 *****************************************************************************/
#if defined(LL_PROFILE_UPDATES)
void LLAMA_DV::updates_start(){
    scoped_lock<shared_mutex_t> xlock(m_lock_checkpoint);
    g_llama_total_time.start();
}

void LLAMA_DV::updates_stop(){
    scoped_lock<shared_mutex_t> xlock(m_lock_checkpoint);
    g_llama_total_time.stop();
}
#endif

/*****************************************************************************
 *                                                                           *
 *  Dump                                                                     *
 *                                                                           *
 *****************************************************************************/

void LLAMA_DV::dump_ostream(std::ostream& out) const {
    shared_lock<shared_mutex_t> cplock(m_lock_checkpoint); // forbid any checkpoint now
    dump_impl(out, *(m_db->graph()));
}

void LLAMA_DV::dump_ostream(std::ostream& out, int level) const {
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

void LLAMA_DV::dump_snapshot(ll_mlcsr_ro_graph& graph) const{
    dump_impl(std::cout, graph);
}

template<typename T>
void LLAMA_DV::dump_impl(std::ostream& out, T& graph) const{
    out << "[LLAMA] num vertices (global): " << num_vertices() << ", num edges (global): " << num_edges() << ", "
            "num levels (#snapshots): " << graph.num_levels() << endl;

    for(node_t node_id = 0; node_id < graph.max_nodes(); node_id++){
        // graph.node_exists is quite inconsistent. Depending whether the node is in the write store or in the read store, the behaviour is different:
        // in the read store, the node exists if has at least one incoming or outgoing edge
        // in the write store, the node exists if it has just been inserted with #insert, regardless of the number of edges it has
        if(!(graph.node_exists(node_id))) continue;

        out << "[" << node_id << "] " << get_write_store_outdegree(node_id) << " outgoing edges: ";

        ll_edge_iterator iterator;
        graph.out_iter_begin(iterator, node_id);
        bool first = true;
        for (edge_t e = graph.out_iter_next(iterator); e != LL_NIL_EDGE; e = graph.out_iter_next(iterator)) {
            node_t n = LL_ITER_OUT_NEXT_NODE(graph, iterator, e);

            if(first) first = false; else out << ", ";
            if(e >= 0){
                out << "<" << n << ", e_id: " << e << ", weight: " << get_out_edge_weight(graph, e) << ">";
            } else {
                out << n;
            }
        }

        if(!m_is_directed){
            graph.in_iter_begin_fast(iterator, node_id);
            for (edge_t e = graph.in_iter_next_fast(iterator); e != LL_NIL_EDGE; e = graph.in_iter_next_fast(iterator)) {
                node_t n = LL_ITER_OUT_NEXT_NODE(graph, iterator, e);

                if(first) first = false; else out << ", ";
                if(e >= 0){
                    out << "<" << n << ", e_id: " << e << ", e_weight: " << get_in_edge_weight(graph, e) << ">";
                } else {
                    out << n;
                }
            }
        }

        out << "\n";
    }
}

// explicitly instantiate the method for both ll_mlcsr_ro_graph and ll_writable_graph
template
void LLAMA_DV::dump_impl<ll_mlcsr_ro_graph>(std::ostream& out, ll_mlcsr_ro_graph& graph) const;
template
void LLAMA_DV::dump_impl<ll_writable_graph>(std::ostream& out, ll_writable_graph& graph) const;

/*****************************************************************************
 *                                                                           *
 *  BFS                                                                      *
 *                                                                           *
 *****************************************************************************/

// llama's implementation
namespace { // anonymous

template<class Graph, typename level_t, bool use_multithread, bool has_navigator,
        bool use_reverse_edge, bool save_child>
class ll_bfs_template {

  public:
    ll_bfs_template(Graph& _G, bool is_graph_undirected) :
        is_undirected(is_graph_undirected), G(_G) {
        visited_bitmap = NULL; // bitmap
        visited_level = NULL;
        thread_local_next_level = NULL;
        down_edge_array = NULL;
        down_edge_set = NULL;
        down_edge_array_w = NULL;
        if (save_child) {
            down_edge_set = new std::unordered_set<edge_t>();
        }
    }

    virtual ~ll_bfs_template() {
        delete [] visited_bitmap;
        delete [] visited_level;
        delete [] thread_local_next_level;
        delete down_edge_set;

                if (down_edge_array != NULL) {
#ifndef FORCE_L0
                        for (size_t i = 0; i < G.num_levels(); i++) delete[] down_edge_array[i];
#endif
                        delete[] down_edge_array;
                }
    }

    void prepare(node_t root_node) {
                // TODO Is this correct? Do we need to poll a some sort of a runtime?
                prepare(root_node, omp_get_max_threads());
        }

    void prepare(node_t root_node, int max_num_thread) {
        int num_thread;
        if (use_multithread) {
            num_thread = max_num_thread;
        } else {
            num_thread = 1;
        }
                max_threads = num_thread;

        is_finished = false;
        curr_level = 0;
        root = root_node;
        state = ST_SMALL;
        assert(root != LL_NIL_NODE);
        if (save_child) {
            if (down_edge_set == NULL)
                down_edge_set = new std::unordered_set<edge_t>();
        }

        global_vector.clear();
        level_queue_begin.clear();
        level_count.clear();
        // create local queues
        if (thread_local_next_level == NULL) {
            thread_local_next_level = new std::vector<node_t>[num_thread];
            for (int i = 0; i < num_thread; i++)
                thread_local_next_level[i].reserve(THRESHOLD2);
        } else {
            for (int i = 0; i < num_thread; i++)
                thread_local_next_level[i].clear();
        }
    }

    void do_bfs_forward(utility::TimeoutService& time_budget) {
        //---------------------------------
        // prepare root node
        //---------------------------------
        curr_level = 0;
        curr_count = 0;
        next_count = 0;

        small_visited[root] = curr_level;
        curr_count++;
        global_vector.push_back(root);
        global_curr_level_begin = 0;
        global_next_level_begin = curr_count;

        level_count.push_back(curr_count);
        level_queue_begin.push_back(global_curr_level_begin);

        bool is_done = false;
        while (!is_done && !time_budget.is_timeout()) {
            switch (state) {
                case ST_SMALL: {
                    for (node_t i = 0; i < curr_count; i++) {
                        node_t t = global_vector[global_curr_level_begin + i];
//                        iterate_neighbor_small(t);
                        visit(t, &ll_bfs_template::visit_neighbor_small);
                        visit_fw(t);            // visit after iteration. in that way, one can check  down-neighbors quite easily
                    }
                    break;
                }
                case ST_QUE: {
                    if (use_multithread)  // do it in parallel
                    {
                        int num_threads = std::min((node_t) max_threads, curr_count/128+1);
                        #pragma omp parallel num_threads(num_threads)
                        {
                            int tid = omp_get_thread_num();
                            #pragma omp for nowait
                            for (node_t i = 0; i < curr_count; i++) {
                                node_t t = global_vector[global_curr_level_begin + i];
//                                iterate_neighbor_que(t, tid);
                                visit(t, &ll_bfs_template::visit_neighbor_que, tid);
                                visit_fw(t);
                            }
                            finish_thread_que(tid);
                        }
                    }
                    else { // do it in sequential
                            int tid = 0;
                            for (node_t i = 0; i < curr_count; i++) {
                                //node_t t = global_curr_level[i];
                                node_t t = global_vector[global_curr_level_begin + i];
//                                iterate_neighbor_que(t, tid);
                                visit(t, &ll_bfs_template::visit_neighbor_que, tid);
                                visit_fw(t);
                            }
                            finish_thread_que(tid);
                    }
                    break;
                }
                case ST_Q2R: {
                    if (use_multithread) {  // do it in parallel
                        int num_threads = std::min((node_t) max_threads, curr_count/128+1);
                        #pragma omp parallel num_threads(num_threads)
                        {
                            node_t local_cnt = 0;
                            #pragma omp for nowait
                            for (node_t i = 0; i < curr_count; i++) {
                                node_t t = global_vector[global_curr_level_begin + i];
//                                iterate_neighbor_rd(t, local_cnt);
                                visit(t, &ll_bfs_template::visit_neighbor_rd, local_cnt);
                                visit_fw(t);
                            }
                            finish_thread_rd(local_cnt);
                        }
                    } else { // do it sequentially
                            node_t local_cnt = 0;
                            for (node_t i = 0; i < curr_count; i++) {
                                //node_t t = global_curr_level[i];
                                node_t t = global_vector[global_curr_level_begin + i];
//                                iterate_neighbor_rd(t, local_cnt);
                                visit(t, &ll_bfs_template::visit_neighbor_rd, local_cnt);
                                visit_fw(t);
                            }
                            finish_thread_rd(local_cnt);
                    }
                    break;
                }

                case ST_RD: {
                    if (use_multithread) { // do it in parallel
                        #pragma omp parallel
                        {
                            node_t local_cnt = 0;
                            #pragma omp for nowait schedule(dynamic,128)
                            for (node_t t = 0; t < G.max_nodes(); t++) {
                                if (visited_level[t] == curr_level) {
//                                    iterate_neighbor_rd(t, local_cnt);
                                    visit(t, &ll_bfs_template::visit_neighbor_rd, local_cnt);
                                    visit_fw(t);
                                }
                            }
                            finish_thread_rd(local_cnt);
                        }
                    } else { // do it in sequential
                            node_t local_cnt = 0;
                            for (node_t t = 0; t < G.max_nodes(); t++) {
                                if (visited_level[t] == curr_level) {
//                                    iterate_neighbor_rd(t, local_cnt);
                                    visit(t, &ll_bfs_template::visit_neighbor_rd, local_cnt);
                                    visit_fw(t);
                                }
                            }
                            finish_thread_rd(local_cnt);
                    }
                    break;
                }
                case ST_R2Q: {
                    if (use_multithread) { // do it in parallel
                        #pragma omp parallel
                        {
                            int tid = omp_get_thread_num();
                            #pragma omp for nowait schedule(dynamic,128)
                            for (node_t t = 0; t < G.max_nodes(); t++) {
                                if (visited_level[t] == curr_level) {
//                                    iterate_neighbor_que(t, tid);
                                    visit(t, &ll_bfs_template::visit_neighbor_que, tid);
                                    visit_fw(t);
                                }
                            }
                            finish_thread_que(tid);
                        }
                    } else {
                            int tid = 0;
                            for (node_t t = 0; t < G.max_nodes(); t++) {
                                if (visited_level[t] == curr_level) {
//                                    iterate_neighbor_que(t, tid);
                                    visit(t, &ll_bfs_template::visit_neighbor_que, tid);
                                    visit_fw(t);
                                }
                            }
                            finish_thread_que(tid);
                    }
                    break;
                }
            } // end of switch

            do_end_of_level_fw();
            is_done = get_next_state();

        } // end of while
    }

    void do_bfs_reverse() {
        // This function should be called only after do_bfs_foward has finished.
        // assumption: small-world graph
        level_t& level = curr_level;
        while (true) {
            node_t count = level_count[level];
            //node_t* queue_ptr = level_start_ptr[level];
            node_t* queue_ptr;
            node_t begin_idx = level_queue_begin[level];
            if (begin_idx == -1) {
                queue_ptr = NULL;
            } else {
                queue_ptr = & (global_vector[begin_idx]);
            }

            if (queue_ptr == NULL) {
#pragma omp parallel if (use_multithread)
                {
#pragma omp for nowait schedule(dynamic,128)
                    for (node_t i = 0; i < G.max_nodes(); i++) {
                        if (visited_level[i] != curr_level) continue;
                        visit_rv(i, level);
                    }
                }
            } else {
                                int num_threads = std::min((node_t) max_threads, curr_count/128+1);
#pragma omp parallel num_threads(num_threads) if (use_multithread)
                {
#pragma omp for nowait
                    for (node_t i = 0; i < count; i++) {
                        node_t u = queue_ptr[i];
                        visit_rv(u, level);
                    }
                }
            }

            do_end_of_level_rv();
            if (level == 0) break;
            level--;
        }
    }

    bool is_down_edge(edge_t idx) {
        if (state == ST_SMALL)
            return (down_edge_set->find(idx) != down_edge_set->end());
        else {
#ifdef FORCE_L0
            return down_edge_array[idx];
#else
                        size_t level = LL_EDGE_LEVEL(idx);
                        if (level == LL_WRITABLE_LEVEL) {
                                return down_edge_array_w[LL_EDGE_GET_WRITABLE(idx)->we_numerical_id];
                        }
                        return down_edge_array[level][LL_EDGE_INDEX(idx)];
#endif
                }
    }

  protected:
    virtual void visit_fw(node_t t)=0;
    virtual void visit_rv(node_t t)=0;
    virtual bool check_navigator(node_t t, edge_t nx)=0;
    virtual void do_end_of_level_fw() {
    }
    virtual void do_end_of_level_rv() {
    }

    node_t get_root() {
        return root;
    }

public:
    level_t get_level(node_t t) {
        // GCC expansion
        if (__builtin_expect((state == ST_SMALL), 0)) {
            if (small_visited.find(t) == small_visited.end())
                return __INVALID_LEVEL;
            else
                return small_visited[t];
        } else {
            return visited_level[t];
        }
    }

private:
    level_t get_curr_level() {
        return curr_level;
    }


  private:
    bool get_next_state() {
        //const char* state_name[5] = {"SMALL","QUEUE","Q2R","RD","R2Q"};

        if (next_count == 0) return true;  // BFS is finished

        int next_state = state;
        switch (state) {
            case ST_SMALL:
                if (next_count >= THRESHOLD1) {
                    prepare_que();
                    next_state = ST_QUE;
                }
                break;
            case ST_QUE:
                if ((next_count >= THRESHOLD2) && (next_count >= curr_count*5)) {
                    prepare_read();
                    next_state = ST_Q2R;
                }
                break;
            case ST_Q2R:
                next_state = ST_RD;
                break;
            case ST_RD:
                if (next_count <= (2 * curr_count)) {
                    next_state = ST_R2Q;
                }
                break;
            case ST_R2Q:
                next_state = ST_QUE;
                break;
        }

        finish_level(state);
        state = next_state;

        return false;
    }

    void finish_level(int state) {
        if ((state == ST_RD) || (state == ST_Q2R)) {
            // output queue is not valid
        } else { // move output queue
            //node_t* temp = &(global_next_level[next_count]);
            //global_curr_level = global_next_level;
            //global_next_level = temp;
            global_curr_level_begin = global_next_level_begin;
            global_next_level_begin = global_next_level_begin + next_count;
        }

        curr_count = next_count;
        next_count = 0;
        curr_level++;

        // save 'new current' level status
        level_count.push_back(curr_count);
        if ((state == ST_RD) || (state == ST_Q2R)) {
            //level_start_ptr.push_back(NULL);
            level_queue_begin.push_back(-1);
        } else {
            //level_start_ptr.push_back(global_curr_level);
            level_queue_begin.push_back(global_curr_level_begin);
        }
    }

        void iter_begin(ll_edge_iterator& iter, node_t v) {
        if (use_reverse_edge) {
                        G.in_iter_begin_fast(iter, v);
        } else {
                        G.out_iter_begin(iter, v);
        }
    }

        edge_t iter_next(ll_edge_iterator& iter) {
        if (use_reverse_edge) {
            return G.in_iter_next_fast(iter);
        } else {
            return G.out_iter_next(iter);
        }
    }

    node_t get_node(ll_edge_iterator& iter) {
                return iter.last_node;
    }

//    void iterate_neighbor_small(node_t t) {
//                ll_edge_iterator iter; iter_begin(iter, t);
//                for (edge_t nx = iter_next(iter); nx != LL_NIL_EDGE; nx = iter_next(iter)) {
//            node_t u = get_node(iter);
//
//            // check visited
//            if (small_visited.find(u) == small_visited.end()) {
//                if (has_navigator) {
//                    if (check_navigator(u, nx) == false) continue;
//                }
//
//                if (save_child) {
//                    save_down_edge_small(nx);
//                }
//
//                small_visited[u] = curr_level + 1;
//                //global_next_level[next_count++] = u;
//                global_vector.push_back(u);
//                next_count++;
//            }
//            else if (save_child) {
//                if (has_navigator) {
//                    if (check_navigator(u, nx) == false) continue;
//                }
//
//                if (small_visited[u] == (curr_level+1)){
//                    save_down_edge_small(nx);
//                }
//            }
//        }
//    }

    void visit_neighbor_small(node_t node_src, edge_t edge_id, node_t node_dst) {
        // check visited
        if (small_visited.find(node_dst) == small_visited.end()) {
            if (has_navigator) {
                if (check_navigator(node_dst, edge_id) == false) return;
            }

            if (save_child) {
                save_down_edge_small(edge_id);
            }

            small_visited[node_dst] = curr_level + 1;
            //global_next_level[next_count++] = u;
            global_vector.push_back(node_dst);
            next_count++;
        }
        else if (save_child) {
            if (has_navigator) {
                if (check_navigator(node_dst, edge_id) == false) return;
            }

            if (small_visited[node_dst] == (curr_level+1)){
                save_down_edge_small(edge_id);
            }
        }
    }

    // should be used only when save_child is enabled
    void save_down_edge_small(edge_t idx) {
        down_edge_set->insert(idx);
    }

    void save_down_edge_large(edge_t idx) {
#ifdef FORCE_L0
        down_edge_array[idx] = 1;
#else
                size_t level = LL_EDGE_LEVEL(idx);
                if (level == LL_WRITABLE_LEVEL) {
                        down_edge_array_w[LL_EDGE_GET_WRITABLE(idx)->we_numerical_id] = 1;
                }
                down_edge_array[LL_EDGE_LEVEL(idx)][LL_EDGE_INDEX(idx)] = 1;
#endif
        }

    void prepare_que() {

        global_vector.reserve(G.max_nodes());

        // create bitmap and edges
        if (visited_bitmap == NULL) {
            visited_bitmap = new unsigned char[(G.max_nodes() + 7) / 8];
            visited_level = new level_t[G.max_nodes()];
        }
        if (save_child) {
            if (down_edge_array == NULL) {
#ifdef FORCE_L0
                down_edge_array = new unsigned char [G.max_edges(0)];
#else
                down_edge_array = new unsigned char* [G.num_levels()];
                                for (size_t i = 0; i < G.num_levels(); i++)
                                        down_edge_array[i] = new unsigned char [G.max_edges(i)];
                                // Note: This makes sense only if the current graph is writable,
                                // but fortunatelly it is never accessed unless we are on the
                                // writable level
                                down_edge_array_w = down_edge_array[G.num_levels() - 1];
#endif
                        }
        }

        if (use_multithread) {
                        #pragma omp parallel
            {
                                #pragma omp for nowait
                for (node_t i = 0; i < (G.max_nodes() + 7) / 8; i++)
                    visited_bitmap[i] = 0;

                                #pragma omp for nowait
                for (node_t i = 0; i < G.max_nodes(); i++)
                    visited_level[i] = __INVALID_LEVEL;

                if (save_child) {
#ifdef FORCE_L0
                                        #pragma omp for nowait
                                        for (edge_t i = 0; i < G.max_edges(0); i++)
                        down_edge_array[i] = 0;
#else
                                        #pragma omp for nowait
                                        for (size_t i = 0; i < G.num_levels(); i++)
                        memset(down_edge_array[i], 0, sizeof(unsigned char) * G.max_edges(i));
#endif
                }
            }
        } else {
            for (node_t i = 0; i < (G.max_nodes() + 7) / 8; i++)
                visited_bitmap[i] = 0;
            for (node_t i = 0; i < G.max_nodes(); i++)
                visited_level[i] = __INVALID_LEVEL;
            if (save_child) {
#ifdef FORCE_L0
                                for (edge_t i = 0; i < G.max_edges(0); i++)
                                        down_edge_array[i] = 0;
#else
                                for (size_t i = 0; i < G.num_levels(); i++)
                                        memset(down_edge_array[i], 0, sizeof(unsigned char) * G.max_edges(i));
#endif
            }
        }

        //typename std::unordered_map<node_t, level_t>::iterator II;
        typename std::map<node_t, level_t>::iterator II;
        for (II = small_visited.begin(); II != small_visited.end(); II++) {
            node_t u = II->first;
            level_t lev = II->second;
            _ll_set_bit(visited_bitmap, u);
            visited_level[u] = lev;
        }

        if (save_child) {
            typename std::unordered_set<edge_t>::iterator J;
            for (J = down_edge_set->begin(); J != down_edge_set->end(); J++) {
                                edge_t idx = *J;
#ifdef FORCE_L0
                                down_edge_array[idx] = 1;
#else
                                size_t level = LL_EDGE_LEVEL(idx);
                                if (level == LL_WRITABLE_LEVEL) {
                                        down_edge_array_w[LL_EDGE_GET_WRITABLE(idx)->we_numerical_id] = 1;
                                }
                                down_edge_array[level][LL_EDGE_INDEX(idx)] = 1;
#endif
            }
        }
    }

//    void iterate_neighbor_que(node_t t, int tid) {
//                ll_edge_iterator iter; iter_begin(iter, t);
//                for (edge_t nx = iter_next(iter); nx != LL_NIL_EDGE; nx = iter_next(iter)) {
//            node_t u = get_node(iter);
//                        assert(u >= 0 && u < G.max_nodes());
//
//            // check visited bitmap
//            // test & test& set
//            if (_ll_get_bit(visited_bitmap, u) == 0) {
//                if (has_navigator) {
//                    if (check_navigator(u, nx) == false) continue;
//                }
//
//                bool re_check_result;
//                if (use_multithread) {
//                    re_check_result = _ll_set_bit_atomic(visited_bitmap, u);
//                } else {
//                    re_check_result = true;
//                    _ll_set_bit(visited_bitmap, u);
//                }
//
//                if (save_child) {
//                    save_down_edge_large(nx);
//                }
//
//                if (re_check_result) {
//                    // add to local q
//                    thread_local_next_level[tid].push_back(u);
//                    visited_level[u] = (curr_level + 1);
//                }
//            }
//            else if (save_child) {
//                if (has_navigator) {
//                    if (check_navigator(u, nx) == false) continue;
//                }
//                if (visited_level[u] == (curr_level +1)) {
//                    save_down_edge_large(nx);
//                }
//            }
//        }
//    }

    void visit_neighbor_que(node_t node_src, edge_t edge_id, node_t node_dst, int tid) {
        assert(node_dst >= 0 && node_dst < G.max_nodes());

        // check visited bitmap
        // test & test& set
        if (_ll_get_bit(visited_bitmap, node_dst) == 0) {
            if (has_navigator) {
                if (check_navigator(node_dst, edge_id) == false) return;
            }

            bool re_check_result;
            if (use_multithread) {
                re_check_result = _ll_set_bit_atomic(visited_bitmap, node_dst);
            } else {
                re_check_result = true;
                _ll_set_bit(visited_bitmap, node_dst);
            }

            if (save_child) {
                save_down_edge_large(edge_id);
            }

            if (re_check_result) {
                // add to local q
                thread_local_next_level[tid].push_back(node_dst);
                visited_level[node_dst] = (curr_level + 1);
            }
        }
        else if (save_child) {
            if (has_navigator) {
                if (check_navigator(node_dst, edge_id) == false) return;
            }
            if (visited_level[node_dst] == (curr_level +1)) {
                save_down_edge_large(edge_id);
            }
        }
    }

    void finish_thread_que(int tid) {
        node_t local_cnt = thread_local_next_level[tid].size();
        //copy curr_cnt to next_cnt
        if (local_cnt > 0) {
            node_t old_idx = __sync_fetch_and_add(&next_count, local_cnt);
            // copy to global vector
            memcpy(&(global_vector[global_next_level_begin + old_idx]),
                   &(thread_local_next_level[tid][0]),
                   local_cnt * sizeof(node_t));
        }
        thread_local_next_level[tid].clear();
    }

    void prepare_read() {
        // nothing to do
    }

//    void iterate_neighbor_rd(node_t t, node_t& local_cnt) {
//                ll_edge_iterator iter; iter_begin(iter, t);
//                for (edge_t nx = iter_next(iter); nx != LL_NIL_EDGE; nx = iter_next(iter)) {
//            node_t u = get_node(iter);
//
//            // check visited bitmap
//            // test & test& set
//            if (_ll_get_bit(visited_bitmap, u) == 0) {
//                if (has_navigator) {
//                    if (check_navigator(u, nx) == false) continue;
//                }
//
//                bool re_check_result;
//                if (use_multithread) {
//                    re_check_result = _ll_set_bit_atomic(visited_bitmap, u);
//                } else {
//                    re_check_result = true;
//                    _ll_set_bit(visited_bitmap, u);
//                }
//
//                if (save_child) {
//                    save_down_edge_large(nx);
//                }
//
//                if (re_check_result) {
//                    // add to local q
//                    visited_level[u] = curr_level + 1;
//                    local_cnt++;
//                }
//            }
//            else if (save_child) {
//                if (has_navigator) {
//                    if (check_navigator(u, nx) == false) continue;
//                }
//                if (visited_level[u] == (curr_level +1)) {
//                    save_down_edge_large(nx);
//                }
//            }
//        }
//    }

    void visit_neighbor_rd(node_t node_src, edge_t edge_id, node_t node_dst, node_t& local_cnt) {
        // check visited bitmap
        // test & test& set
        if (_ll_get_bit(visited_bitmap, node_dst) == 0) {
            if (has_navigator) {
                if (check_navigator(node_dst, edge_id) == false) return;
            }

            bool re_check_result;
            if (use_multithread) {
                re_check_result = _ll_set_bit_atomic(visited_bitmap, node_dst);
            } else {
                re_check_result = true;
                _ll_set_bit(visited_bitmap, node_dst);
            }

            if (save_child) {
                save_down_edge_large(edge_id);
            }

            if (re_check_result) {
                // add to local q
                visited_level[node_dst] = curr_level + 1;
                local_cnt++;
            }
        }
        else if (save_child) {
            if (has_navigator) {
                if (check_navigator(node_dst, edge_id) == false) return;
            }
            if (visited_level[node_dst] == (curr_level +1)) {
                save_down_edge_large(edge_id);
            }
        }
    }

    template<typename Function, typename ...Args>
    void visit(node_t node_src, Function lambda, Args&&... args){
        ll_edge_iterator iterator;
        G.out_iter_begin(iterator, node_src);
        for (edge_t edge_id = G.out_iter_next(iterator); edge_id != LL_NIL_EDGE; edge_id = G.out_iter_next(iterator)) {
            node_t node_dst = get_node(iterator);
            std::invoke(lambda, this, node_src, edge_id, node_dst, std::forward<Args>(args)...);
        }

        // if the graph is undirected, then the actual set of  edges is comprised by both the stored outgoing and incoming edges
        if(is_undirected){
            G.in_iter_begin_fast(iterator, node_src);
            for (edge_t edge_id = G.in_iter_next_fast(iterator); edge_id != LL_NIL_EDGE; edge_id = G.in_iter_next_fast(iterator)) {
                node_t node_dst = get_node(iterator);
                std::invoke(lambda, this, node_src, edge_id, node_dst, std::forward<Args>(args)...);
            }
        }
    }

    void finish_thread_rd(node_t local_cnt) {
        __sync_fetch_and_add(&next_count, local_cnt);
    }

public:
    //-----------------------------------------------------
    //-----------------------------------------------------
    static const int ST_SMALL = 0;
    static const int ST_QUE = 1;
    static const int ST_Q2R = 2;
    static const int ST_RD = 3;
    static const int ST_R2Q = 4;
    static const int THRESHOLD1 = 128;  // single threaded
    static const int THRESHOLD2 = 1024; // move to RD-based

    // not -1.
    //(why? because curr_level-1 might be -1, when curr_level = 0)
    static const level_t __INVALID_LEVEL = -2;

private:
    int state;

    // whether the graph is undirected or not. If the graph is directed, at each step, only visit the outgoing edges.
    // Otherwise, iterator over both the incoming the and outgoing edges when visiting a node.
    const bool is_undirected;

    unsigned char* visited_bitmap; // bitmap
    level_t* visited_level; // assumption: small_world graph
    bool is_finished;
    level_t curr_level;
    node_t root;
    Graph& G;
    node_t curr_count;
    node_t next_count;

    //std::unordered_map<node_t, level_t> small_visited;
    std::map<node_t, level_t> small_visited;
    std::unordered_set<edge_t>* down_edge_set;
        unsigned char* down_edge_array_w;
#ifdef FORCE_L0
    unsigned char* down_edge_array;
#else
    unsigned char** down_edge_array;
#endif

    //node_t* global_next_level;
    //node_t* global_curr_level;
    //node_t* global_queue;
    std::vector<node_t> global_vector;
    node_t global_curr_level_begin;
    node_t global_next_level_begin;

    //std::vector<node_t*> level_start_ptr;
    std::vector<node_t> level_queue_begin;
    std::vector<node_t> level_count;

    std::vector<node_t>* thread_local_next_level;

        int max_threads;
};


template <class Graph>
class bfs_bfs : public ll_bfs_template
    <Graph, int16_t, true, false, false, false>
{
public:
    bfs_bfs(Graph& _G, node_t& _root, bool is_graph_undirected)
    : ll_bfs_template<Graph, int16_t, true, false, false, false>(_G, is_graph_undirected),
    G(_G), root(_root){}

private:  // list of varaibles
    Graph& G;
    node_t& root;

public:

protected:
    virtual void visit_fw(node_t v){ }

    virtual void visit_rv(node_t v) {}
    virtual bool check_navigator(node_t v, edge_t v_idx) {return true;}


};

} // anonymous namespace


void LLAMA_DV::bfs(uint64_t ext_vertex_id, const char* dump2file) {
    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();

    auto graph = get_snapshot();

    // execute the BFS algorithm
    int64_t vertex_id = ext_vertex_id;
    bfs_bfs<ll_mlcsr_ro_graph> instance{ graph, vertex_id, is_undirected() };
    instance.prepare(vertex_id);
    instance.do_bfs_forward(timeout);

    if(timeout.is_timeout()){
        RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);
    }

    // store the results in the given file
    if(dump2file != nullptr){
        COUT_DEBUG("save the results to: " << dump2file)
        fstream handle(dump2file, ios_base::out);
        if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

        for(node_t node_id = 0; node_id < graph.max_nodes(); node_id++){
            // first, does this node exist (or it's a gap?)
            // this is a bit of a stretch: the impl~ from llama assumes that a node does not exist only if it does not have any incoming or outgoing edges.
            if(!graph.node_exists(node_id)) continue;
            handle << node_id << " ";

            int distance = instance.get_level(node_id);

            if(distance != decltype(instance)::__INVALID_LEVEL){
                handle << distance;
            } else {
                handle << std::numeric_limits<int64_t>::max(); // it should have been -1, but ok
            }

            handle << "\n";

        }
        handle.close();
    }
}

/*****************************************************************************
 *                                                                           *
 *  CDLP                                                                     *
 *                                                                           *
 *****************************************************************************/

void LLAMA_DV::cdlp(uint64_t max_iterations, const char* dump2file){
    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();

    auto graph = get_snapshot();

    // execute the CDLP algortihm
    unique_ptr<uint64_t[]> ptr_labels = cdlp_impl(timeout, graph, max_iterations);
    uint64_t* labels = ptr_labels.get();

    if(timeout.is_timeout()){
        RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);
    }

    // store the results in the given file
    if(dump2file != nullptr){
        COUT_DEBUG("save the results to: " << dump2file)
        fstream handle(dump2file, ios_base::out);
        if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

        for(node_t node_id = 0; node_id < graph.max_nodes(); node_id++){
            // first, does this node exist (or it's a gap?)
            // this is a bit of a stretch: the impl~ from llama assumes that a node does not exist only if it does not have any incoming or outgoing edges.
            if(!graph.node_exists(node_id)) continue;

            handle << node_id << " " << labels[node_id] << "\n";
        }

        handle.close();
    }
}

unique_ptr<uint64_t[]> LLAMA_DV::cdlp_impl(utility::TimeoutService& timer, ll_mlcsr_ro_graph& graph, uint64_t max_iterations){
    unique_ptr<uint64_t[]> ptr_labels0 { new uint64_t[graph.max_nodes()] };
    unique_ptr<uint64_t[]> ptr_labels1 { new uint64_t[graph.max_nodes()] };
    uint64_t* labels0 = ptr_labels0.get(); // current labels
    uint64_t* labels1 = ptr_labels1.get(); // labels for the next iteration

    // initialisation
    #pragma omp parallel for
    for(node_t n = 0; n < graph.max_nodes(); n++){
        labels0[n] = n;
    }

    // algorithm pass
    bool change = true;
    uint64_t current_iteration = 0;
    while(current_iteration < max_iterations && change && !timer.is_timeout()){
        change = false; // reset the flag

        #pragma omp parallel for shared(change)
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

/*****************************************************************************
 *                                                                           *
 *  LCC                                                                      *
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
void llama_execute_lcc_directed(utility::TimeoutService& timer, ll_mlcsr_ro_graph& graph, double* scores){

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
void llama_execute_lcc_undirected(utility::TimeoutService& timer, ll_mlcsr_ro_graph& graph, double* scores){

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

void LLAMA_DV::lcc(const char* dump2file){
    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();

    auto graph = get_snapshot();

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

    // store the results in the given file
    if(dump2file != nullptr){
        COUT_DEBUG("save the results to: " << dump2file)
        fstream handle(dump2file, ios_base::out);
        if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

        for(node_t node_id = 0; node_id < graph.max_nodes(); node_id++){
            // first, does this node exist (or it's a gap?)
            // this is a bit of a stretch: the impl~ from llama assumes that a node does not exist only if it does not have any incoming or outgoing edges.
            if(!graph.node_exists(node_id)) continue;

            handle << node_id << " " << scores[node_id] << "\n";
        }

        handle.close();
    }
}

/*****************************************************************************
 *                                                                           *
 *  Pagerank                                                                 *
 *                                                                           *
 *****************************************************************************/

void LLAMA_DV::pagerank(uint64_t num_iterations, double damping_factor, const char* dump2file){
    utility::TimeoutService timeout_srv { m_timeout };
    Timer timer; timer.start();

    auto graph = get_snapshot();
    auto current_num_vertices = num_vertices();

    // execute the Pagerank algorithm
    unique_ptr<double[]> ptr_rank { new double[graph.max_nodes()] };
    double* rank = ptr_rank.get();
    pagerank_impl(timeout_srv, graph, current_num_vertices, num_iterations, damping_factor, /* output */ rank);

    if(timeout_srv.is_timeout()){
        RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);
    }

    // store the results in the given file
    if(dump2file != nullptr){
        COUT_DEBUG("save the results to: " << dump2file)
        fstream handle(dump2file, ios_base::out);
        if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

        for(node_t node_id = 0; node_id < graph.max_nodes(); node_id++){
            // first, does this node exist (or it's a gap?)
            // this is a bit of a stretch: the impl~ from llama assumes that a node does not exist only if it does not have any incoming or outgoing edges.
            if(!graph.node_exists(node_id)) continue;

            handle << node_id << " " << rank[node_id] << "\n";
        }

        handle.close();
    }
}

// Implementation derived from llama/benchmark/benchmarks/pagerank.h, class ll_b_pagerank_pull_ext
void LLAMA_DV::pagerank_impl(utility::TimeoutService& timer, ll_mlcsr_ro_graph& graph, uint64_t num_vertices, uint64_t num_iterations, double d, /* output array, already allocated */ double* G_pg_rank){
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

/*****************************************************************************
 *                                                                           *
 *  SSSP                                                                     *
 *                                                                           *
 *****************************************************************************/
// the following source code is based on the class `ll_b_sssp_weighted', located
// in llama/benchmark/benchmarks/sssp.h
// @param G_dist is the output array with the distances, it must be already allocated (but not initialised) and its length equal at least to graph.max_nodes()
static void llama_execute_sssp(utility::TimeoutService& timer, ll_mlcsr_ro_graph& graph, node_t root, const char* weights_property_name, bool is_graph_undirected, double* __restrict G_dist){
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

void LLAMA_DV::sssp(uint64_t vertex_id, const char* dump2file){
    utility::TimeoutService timeout_srv { m_timeout };
    Timer timer; timer.start();

    auto graph = get_snapshot();

    // execute the SSSP algorithm
    unique_ptr<double[]> ptr_distances { new double[graph.max_nodes()] };
    double* distances = ptr_distances.get();
    llama_execute_sssp(timeout_srv, graph, vertex_id, g_llama_property_weights, is_undirected(), /* output */ distances);

    if(timeout_srv.is_timeout()){
        RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);
    }

    // store the results in the given file
    if(dump2file != nullptr){
        COUT_DEBUG("save the results to: " << dump2file)
        fstream handle(dump2file, ios_base::out);
        if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

        for(node_t node_id = 0; node_id < graph.max_nodes(); node_id++){
            // first, does this node exist (or it's a gap?)
            // this is a bit of a stretch: the impl~ from llama assumes that a node does not exist only if it does not have any incoming or outgoing edges.
            if(!graph.node_exists(node_id)) continue;

            handle << node_id << " " << distances[node_id] << "\n";
        }

        handle.close();
    }
}

/*****************************************************************************
 *                                                                           *
 *  WCC                                                                      *
 *                                                                           *
 *****************************************************************************/

namespace { // anonymous

template<class Graph, bool has_pre_visit, bool has_post_visit, bool has_navigator>
class ll_dfs_template {
protected:
    utility::TimeoutService& timeout_srv;
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
    ll_dfs_template(utility::TimeoutService& timeout_srv, Graph& _G) : timeout_srv(timeout_srv), G(_G) {
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
    WeaklyConnectedComponents(utility::TimeoutService& timeout_srv, Graph& _G, node_t*& _G_SCC, node_t _n)
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
void llama_execute_wcc(utility::TimeoutService& timeout_srv, ll_mlcsr_ro_graph& graph, /* output array, already allocated */ node_t* G_SCC){

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

void LLAMA_DV::wcc(const char* dump2file){
    utility::TimeoutService timeout_srv { m_timeout };
    Timer timer; timer.start();

    auto graph = get_snapshot();

    // execute the WCC algorithm
    unique_ptr<node_t[]> ptr_components { new node_t[graph.max_nodes()] };
    node_t* components = ptr_components.get();
    llama_execute_wcc(timeout_srv, graph, /* output */ components);

    if(timeout_srv.is_timeout()){
        RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);
    }

    // store the results in the given file
    if(dump2file != nullptr){
        COUT_DEBUG("save the results to: " << dump2file)
        fstream handle(dump2file, ios_base::out);
        if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

        for(node_t node_id = 0; node_id < graph.max_nodes(); node_id++){
            // first, does this node exist (or it's a gap?)
            // this is a bit of a stretch: the impl~ from llama assumes that a node does not exist only if it does not have any incoming or outgoing edges.
            if(!graph.node_exists(node_id)) continue;

            handle << node_id << " " << components[node_id] << "\n";
        }

        handle.close();
    }
}

} // namespace

