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
#include "teseo_driver.hpp"

#include <cassert>
#include <fstream>
#include <iomanip>
#include <mutex>

#include "common/timer.hpp"
#include "third-party/gapbs/gapbs.hpp"
#include "utility/timeout_service.hpp"

#include "teseo/context/global_context.hpp"
#include "teseo/memstore/memstore.hpp"
#include "teseo.hpp"


using namespace common;
using namespace std;
using namespace teseo;

#define TESEO reinterpret_cast<Teseo*>(m_pImpl)

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
#define DEBUG
namespace gfe { extern mutex _log_mutex [[maybe_unused]]; }
#define COUT_DEBUG_FORCE(msg) { std::scoped_lock<std::mutex> lock{::gfe::_log_mutex}; std::cout << "[TeseoDriver::" << __FUNCTION__ << "] " << msg << std::endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif

namespace gfe::library {

/*****************************************************************************
 *                                                                           *
 *  Init                                                                     *
 *                                                                           *
 *****************************************************************************/
TeseoDriver::TeseoDriver(bool is_directed) : m_pImpl(new Teseo()), m_is_directed(is_directed) {
    if(is_directed == true){ throw std::invalid_argument("Only undirected graphs are currently supported by the front-end"); }
}

TeseoDriver::~TeseoDriver(){
    delete TESEO; m_pImpl = nullptr;
}

void TeseoDriver::on_thread_init(int thread_id){
    TESEO->register_thread();
}

void TeseoDriver::on_thread_destroy(int thread_id){
    TESEO->unregister_thread();
}

/*****************************************************************************
 *                                                                           *
 *  Updates & point look ups                                                 *
 *                                                                           *
 *****************************************************************************/

uint64_t TeseoDriver::num_edges() const {
    auto tx = TESEO->start_transaction();
    return tx.num_edges();
}

uint64_t TeseoDriver::num_vertices() const {
    auto tx = TESEO->start_transaction();
    return tx.num_vertices();
}

bool TeseoDriver::has_vertex(uint64_t vertex_id) const {
    auto tx = TESEO->start_transaction();
    return tx.has_vertex(vertex_id);
}

double TeseoDriver::get_weight(uint64_t source, uint64_t destination) const {
    auto tx = TESEO->start_transaction();
    try {
        return tx.get_weight(source, destination);
    } catch (LogicalError& e) {
        return numeric_limits<double>::signaling_NaN();
    }
}

bool TeseoDriver::is_directed() const {
    return m_is_directed;
}

void TeseoDriver::set_timeout(uint64_t seconds){
    m_timeout = chrono::seconds{ seconds };
}

bool TeseoDriver::add_vertex(uint64_t vertex_id){
    auto tx = TESEO->start_transaction();
    try {
        tx.insert_vertex(vertex_id);
        tx.commit();
        return true;
    } catch( LogicalError& e ){
        return false;
    } catch( TransactionConflict& e) {
        return false;
    }
}

bool TeseoDriver::remove_vertex(uint64_t vertex_id){
    COUT_DEBUG("remove vertex: " << vertex_id);
    auto tx = TESEO->start_transaction();
    try {
        tx.remove_vertex(vertex_id);
        tx.commit();
        return true;
    } catch( LogicalError& e ){
        return false;
    } catch( TransactionConflict& e) {
        return false;
    }
}

bool TeseoDriver::add_edge(gfe::graph::WeightedEdge e) {
    auto tx = TESEO->start_transaction();
    try {
        tx.insert_edge(e.source(), e.destination(), e.weight());
        tx.commit();
        return true;
    } catch( LogicalError& e ){
        return false;
    } catch( TransactionConflict& e) {
        return false;
    }
}

bool TeseoDriver::remove_edge(gfe::graph::Edge e) {
    auto tx = TESEO->start_transaction();
    try {
        tx.remove_edge(e.source(), e.destination());
        tx.commit();
        return true;
    } catch( LogicalError& e ){
        return false;
    } catch( TransactionConflict& e) {
        return false;
    }
}

/*****************************************************************************
 *                                                                           *
 *  Dump                                                                     *
 *                                                                           *
 *****************************************************************************/

void TeseoDriver::dump_ostream(std::ostream& out) const {
    auto memstore = teseo::context::global_context()->memstore();
    memstore->dump();
}

/*****************************************************************************
 *                                                                           *
 *  BFS                                                                      *
 *                                                                           *
 *****************************************************************************/
void TeseoDriver::bfs(uint64_t source_vertex_id, const char* dump2file){
    assert(0 && "To be implemented");
}

/*****************************************************************************
 *                                                                           *
 *  PageRank                                                                 *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG_PAGERANK
#if defined(DEBUG_PAGERANK)
#define COUT_DEBUG_PAGERANK(msg) COUT_DEBUG(msg)
#else
#define COUT_DEBUG_PAGERANK(msg)
#endif

// Implementation based on the reference PageRank for the GAP Benchmark Suite
// https://github.com/sbeamer/gapbs
// The reference implementation has been written by Scott Beamer
//
// Copyright (c) 2015, The Regents of the University of California (Regents)
// All Rights Reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
// 3. Neither the name of the Regents nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL REGENTS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

/*
GAP Benchmark Suite
Kernel: PageRank (PR)
Author: Scott Beamer

Will return pagerank scores for all vertices once total change < epsilon

This PR implementation uses the traditional iterative approach. This is done
to ease comparisons to other implementations (often use same algorithm), but
it is not necesarily the fastest way to implement it. It does perform the
updates in the pull direction to remove the need for atomics.
*/

static
unique_ptr<double[]> teseo_pagerank(teseo::Teseo* teseo, teseo::Transaction& transaction, uint64_t num_iterations, double damping_factor, utility::TimeoutService& timer) {
    // init
    const uint64_t num_vertices = transaction.num_vertices();
    COUT_DEBUG("num vertices: " << num_vertices);

    auto iterator = transaction.iterator();
    const double init_score = 1.0 / num_vertices;
    const double base_score = (1.0 - damping_factor) / num_vertices;

    unique_ptr<double[]> ptr_scores{ new double[num_vertices]() }; // avoid memory leaks
    double* scores = ptr_scores.get();
    #pragma omp parallel for
    for(uint64_t v = 0; v < num_vertices; v++){
        scores[v] = init_score;
    }
    gapbs::pvector<double> outgoing_contrib(num_vertices, 0.0);

    // pagerank iterations
    for(uint64_t iteration = 0; iteration < num_iterations && !timer.is_timeout(); iteration++){
        double dangling_sum = 0.0;

        // for each node, precompute its contribution to all of its outgoing neighbours and, if it's a sink,
        // add its rank to the `dangling sum' (to be added to all nodes).

        #pragma omp parallel
        {
            teseo->register_thread();

           #pragma omp for reduction(+:dangling_sum)
           for(uint64_t v = 0; v < num_vertices; v++){
               //uint64_t out_degree = view->get_degree_out(v);
               uint64_t out_degree = transaction.degree(v, /* logical */ true);
               if(out_degree == 0){ // this is a sink
                   dangling_sum += scores[v];
               } else {
                   outgoing_contrib[v] = scores[v] / out_degree;
               }
           }

           teseo->unregister_thread();
        }



        dangling_sum /= num_vertices;

        // compute the new score for each node in the graph
        #pragma omp parallel
        {
            teseo->register_thread();

            #pragma omp for schedule(dynamic, 64)
            for(uint64_t v = 0; v < num_vertices; v++){

                double incoming_total = 0;
                iterator.edges(v, /* logical ? */ true, [&incoming_total, &outgoing_contrib](uint64_t destination, double weight){
                   incoming_total += outgoing_contrib[destination];
                   return true;
                });

                // update the score
                scores[v] = base_score + damping_factor * (incoming_total + dangling_sum);
            }


            teseo->unregister_thread();
        }
    }

    return ptr_scores;
}

void TeseoDriver::pagerank(uint64_t num_iterations, double damping_factor, const char* dump2file) {
    // Init
    utility::TimeoutService timeout { m_timeout };
    Timer timer; timer.start();

    // Run the PageRank algorithm
    auto transaction = TESEO->start_transaction(/* read only ? */ true);
    unique_ptr<double[]> ptr_rank = teseo_pagerank(TESEO, transaction, num_iterations, damping_factor, timeout);
    if(timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer);  }

    // store the results in the given file
    if(dump2file != nullptr){
        TESEO->register_thread();

        COUT_DEBUG("save the results to: " << dump2file);
        fstream handle(dump2file, ios_base::out);
        double* __restrict rank = ptr_rank.get();
        if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

        for(uint64_t i = 0, N = transaction.num_vertices(); i < N; i++){
            handle << transaction.vertex_id(i) << " " << rank[i] << "\n";
        }

        handle.close();
    }
}

/*****************************************************************************
 *                                                                           *
 *  WCC                                                                      *
 *                                                                           *
 *****************************************************************************/
void TeseoDriver::wcc(const char* dump2file) {
    assert(0 && "To be implemented");
}

/*****************************************************************************
 *                                                                           *
 *  CDLP                                                                      *
 *                                                                           *
 *****************************************************************************/
void TeseoDriver::cdlp(uint64_t max_iterations, const char* dump2file) {
    assert(0 && "To be implemented");
}

/*****************************************************************************
 *                                                                           *
 *  LCC                                                                      *
 *                                                                           *
 *****************************************************************************/
void TeseoDriver::lcc(const char* dump2file) {
    assert(0 && "To be implemented");
}

/*****************************************************************************
 *                                                                           *
 *  SSSP                                                                     *
 *                                                                           *
 *****************************************************************************/
void TeseoDriver::sssp(uint64_t source_vertex_id, const char* dump2file) {
    assert(0 && "To be implemented");
}

} // namespace
