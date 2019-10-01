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

#include "aging2_master.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <mutex>
#include <thread>
#include <utility>

#include "common/error.hpp"
#include "common/system.hpp"
#include "common/timer.hpp"
#include "experiment/aging2_experiment.hpp"
#include "experiment/aging2_result.hpp"
#include "graph/edge_stream.hpp"
#include "library/interface.hpp"
#include "aging2_worker.hpp"
#include "build_thread.hpp"
#include "configuration.hpp"

using namespace common;
using namespace std;

/*****************************************************************************
 *                                                                           *
 * Debug                                                                     *
 *                                                                           *
 *****************************************************************************/
extern mutex _log_mutex [[maybe_unused]];
//#define DEBUG

#define COUT_DEBUG_FORCE(msg) { scoped_lock<mutex> lock(_log_mutex); cout << "[Aging2Master::" << __FUNCTION__ << "] [" << concurrency::get_thread_id() << "] " << msg << endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif

/*****************************************************************************
 *                                                                           *
 * Init                                                                      *
 *                                                                           *
 *****************************************************************************/
namespace experiment::details {

Aging2Master::Aging2Master(const Aging2Experiment& parameters) :
    m_parameters(parameters),
    m_vertices_final(move(* (parameters.m_stream->vertex_table().get()) ) ),
    m_is_directed(m_parameters.m_library->is_directed()),
    m_results(parameters) {

    // 1024 is a hack to avoid issues with small graphs
    m_reported_times = new uint64_t[ std::max<uint64_t>( ::ceil (static_cast<double>(num_operations_total()) / m_parameters.m_stream->num_edges() ) + 1, 1024 ) ]();

    m_parameters.m_library->on_main_init(m_parameters.m_num_threads + /* this + builder service */ 2);

    init_cumulative_frequencies();
    init_workers();
    m_parameters.m_library->on_thread_init(m_parameters.m_num_threads);
}

Aging2Master::~Aging2Master(){
    for(auto w: m_workers){ delete w; }
    m_workers.clear();
    m_parameters.m_library->on_thread_destroy(m_parameters.m_num_threads);
    m_parameters.m_library->on_main_destroy();

    delete[] m_vertices_freq; m_vertices_freq = nullptr; m_vertices_freq_sz = 0;

    delete[] m_reported_times; m_reported_times = nullptr;

    delete[] m_vertices2remove; m_vertices2remove = nullptr; // avoid potential leaks due to exceptions
}

void Aging2Master::init_cumulative_frequencies(){
    Timer timer; timer.start();
    LOG("[Aging2] Computing the cumulative frequency for the edges in the input graph ... ");

    // compute the frequencies for the vertices
    m_vertices_freq_sz = std::max<uint64_t>( m_vertices_final.size() * parameters().m_ef_vertices, m_vertices_final.size() );
    m_vertices_freq = new RankFrequency[m_vertices_freq_sz]();
    { // restrict the scope
        auto lst_vertices = m_vertices_final.lock_table();
        int64_t pos = 0;
        for(auto& pair: lst_vertices){ m_vertices_freq[ pos++ ] = { pair.first, pair.second }; }
        lst_vertices.unlock();
    }
    std::sort(m_vertices_freq, m_vertices_freq + m_vertices_final.size(), [](const RankFrequency& r1, const RankFrequency& r2){
        return r1.m_frequency > r2.m_frequency; // decreasing order
    });
    assert(m_vertices_freq_sz >= m_vertices_final.size());
    uint64_t vertices_to_insert = m_vertices_freq_sz - m_vertices_final.size();
    if(vertices_to_insert > 0){
        uint64_t vertex_id = 0;

        int64_t pos_tail = m_vertices_freq_sz -1;
        int64_t pos_head = m_vertices_final.size() -1;
        uint64_t remaining_free_spots = vertices_to_insert;
        while(remaining_free_spots > 0 && pos_tail > 0){
            assert(pos_head >= 0);
            if(remaining_free_spots *  m_vertices_freq_sz >= vertices_to_insert * pos_tail){
                // generate the ID of the vertex to insert
                do { vertex_id ++; } while(m_vertices_final.contains(vertex_id));
                // interpolate the frequency w.r.t. the two neighbours
                uint64_t vertex_freq = m_vertices_freq[pos_head].m_frequency;
                if(pos_tail < m_vertices_freq_sz -1){
                    vertex_freq = (vertex_freq + m_vertices_freq[pos_tail +1].m_frequency) /2;
                }
                m_vertices_freq[pos_tail] = {vertex_id, vertex_freq};
                remaining_free_spots--;
//                COUT_DEBUG("vertex[temp][" << pos_tail << "]: " << vertex_id << ", " << vertex_freq);
            } else {
                m_vertices_freq[pos_tail] = m_vertices_freq[pos_head];
                pos_head--;
//                COUT_DEBUG("vertex[final][" << pos_tail << "]: " << m_vertices_freq[pos_tail].m_vertex_id << ", " << m_vertices_freq[pos_tail].m_frequency << ", pos_head: " << pos_head);
            }

            pos_tail--;
        }
    }

    // okay, finally compute the prefix sum
    for(uint64_t i = 1; i < m_vertices_freq_sz; i++){
        m_vertices_freq[i].m_frequency += m_vertices_freq[i -1].m_frequency;
    }

    // build the index
    uint64_t num_entries = std::max<uint64_t>(1, m_vertices_freq_sz / m_vertices_freq_index_leaf_sz);
    m_vertices_freq_index.rebuild(num_entries);
    for(uint64_t i = 0; i < num_entries; i++){
        m_vertices_freq_index.set_separator_key(i, m_vertices_freq[i * m_vertices_freq_index_leaf_sz].m_frequency);
    }

    timer.stop();
    LOG("[Aging2] Cumulative frequencies computed in " << timer);
}

void Aging2Master::init_workers() {
    Timer timer; timer.start();
    LOG("[Aging2] Initialising " << parameters().m_num_threads << " worker threads ... ");

    m_workers.reserve(parameters().m_num_threads);
    for(uint64_t worker_id = 0; worker_id < parameters().m_num_threads; worker_id++){
        m_workers.push_back ( new Aging2Worker(*this, worker_id) );
    }

    LOG("[Aging2] Workers initialised in " << timer);
}

/*****************************************************************************
 *                                                                           *
 * Experiment                                                                *
 *                                                                           *
 *****************************************************************************/
void Aging2Master::do_run_experiment(){
    LOG("[Aging2] Experiment started ...");
    m_last_progress_reported = 0;
    m_last_time_reported = 0; m_time_start = chrono::steady_clock::now();

    // init the build service (the one that creates the new snapshots/deltas)
    BuildThread build_service { parameters().m_library , static_cast<int>(parameters().m_num_threads) +1, parameters().m_build_frequency };

    Timer timer; timer.start();
    for(auto w: m_workers) w->execute_updates();
    for(auto w: m_workers) w->wait();
    build_service.stop();
    m_parameters.m_library->build(); // flush last changes
    timer.stop();
    LOG("[Aging2] Experiment completed!");
    LOG("[Aging2] Updates performed with " << parameters().m_num_threads << " threads in " << timer);
    m_results.m_completion_time = timer.microseconds();
    m_results.m_num_build_invocations = build_service.num_invocations();
}

void Aging2Master::do_remove_artificial_vertices(){
    Timer timer; timer.start();

    auto lst_vertices = m_vertices_present.lock_table();
    m_vertices2remove = new uint64_t[lst_vertices.size()];
    m_results.m_num_artificial_vertices = 0;
    for(auto vertex : lst_vertices){
        if(!m_vertices_final.contains(vertex.first)){
            m_vertices2remove[m_results.m_num_artificial_vertices] = vertex.first;
            m_results.m_num_artificial_vertices++;
        }
    }
    lst_vertices.unlock();

    for(auto w: m_workers) w->remove_artificial_vertices();
    for(auto w: m_workers) w->wait();
    m_parameters.m_library->build();

    delete[] m_vertices2remove; m_vertices2remove = nullptr;

    LOG("[Aging2] Number of extra vertices: " << m_results.m_num_artificial_vertices << ", "
            "expansion factor: " << static_cast<double>(m_results.m_num_artificial_vertices + m_vertices_final.size()) / m_vertices_final.size());
    timer.stop();
    LOG("[Aging2] Artificial vertices removed in " << timer);
}

Aging2Result Aging2Master::execute(){
    do_run_experiment();
    do_remove_artificial_vertices();

    store_results();
    log_num_vtx_edges();

    return m_results;
}

/*****************************************************************************
 *                                                                           *
 * Utility methods                                                           *
 *                                                                           *
 *****************************************************************************/

uint64_t Aging2Master::distr_max_value() const {
    return m_vertices_freq[m_vertices_freq_sz -1].m_frequency;
}

uint64_t Aging2Master::distr_vertex_offset(uint64_t cumulative_frequency) const {
    assert(cumulative_frequency < distr_max_value());
    uint64_t offset = m_vertices_freq_index.find_lt(cumulative_frequency);
    assert(offset < m_vertices_freq_sz && "Starting position outside of bounds");
    while(cumulative_frequency >= m_vertices_freq[offset].m_frequency) offset++;
    assert(offset < m_vertices_freq_sz && "Otherwise cumulative_frequency >= distr_max_value()");
    return offset;
}

uint64_t Aging2Master::num_operations_total() const {
    return parameters().m_mult_ops * parameters().m_stream->num_edges();
}

uint64_t Aging2Master::max_num_edges() const{
    return parameters().m_ef_edges * parameters().m_stream->num_edges();
}

void Aging2Master::store_results(){
    m_results.m_num_vertices_load = m_vertices_final.size();
    m_results.m_num_edges_load = parameters().m_stream->num_edges();
    m_results.m_num_vertices_final_graph = parameters().m_library->num_vertices();
    m_results.m_num_edges_final_graph = parameters().m_library->num_edges();
    m_results.m_num_operations_expected = num_operations_total();
    m_results.m_num_operations_performed = m_num_operations_performed;
    m_results.m_reported_times.reserve(m_last_time_reported);
    for(size_t i = 0, sz = m_last_time_reported; i < sz; i++){
        m_results.m_reported_times.push_back( m_reported_times[i] );
    }
}

void Aging2Master::log_num_vtx_edges(){
    uint64_t num_vertices0 = m_vertices_final.size();
    uint64_t num_edges0 = parameters().m_stream->num_edges();

    scoped_lock<mutex> lock(_log_mutex);
    cout << "[Aging2] Number of stored vertices: " << m_results.m_num_vertices_final_graph << " [match: ";
    if(num_vertices0 == m_results.m_num_vertices_final_graph){ cout << "yes"; } else {
        cout << "no, expected " << num_edges0;
    }
    cout << "], number of stored edges: " << m_results.m_num_edges_final_graph << " [match: ";
    if(num_edges0 == m_results.m_num_edges_final_graph){ cout << "yes"; } else {
        cout << "no, expected " << num_edges0;
    }
    cout << "]" << endl;
}

} // namespace
