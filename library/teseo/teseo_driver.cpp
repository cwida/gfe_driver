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
#include <mutex>


#include "teseo/context/global_context.hpp"
#include "teseo/memstore/memstore.hpp"
#include "teseo.hpp"

using namespace std;
using namespace teseo;

#define TESEO reinterpret_cast<Teseo*>(m_pImpl)

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG
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

} // namespace
