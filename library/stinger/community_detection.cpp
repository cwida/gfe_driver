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
#include "stinger.hpp"

#include <algorithm>
#include <chrono>
#include <cinttypes>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <mutex>
#include <thread>
#include <unordered_map>

#include "common/system.hpp"
#include "stinger_core/stinger.h"
#include "stinger_core/xmalloc.h"
#include "stinger_error.hpp"

using namespace std;

// Macros
#define STINGER reinterpret_cast<struct stinger*>(m_stinger_graph)
#undef CURRENT_ERROR_TYPE
#define CURRENT_ERROR_TYPE ::gfe::library::StingerError
#define TIMER_INIT auto time_start = chrono::steady_clock::now();
#define CHECK_TIMEOUT if(m_timeout > 0 && (chrono::steady_clock::now() - time_start) > chrono::seconds(m_timeout)) { \
        RAISE_EXCEPTION(TimeoutError, "Timeout occurred after: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - time_start).count() << " seconds") };


/******************************************************************************
 *                                                                            *
 *  Debug                                                                     *
 *                                                                            *
 *****************************************************************************/
//#define DEBUG
extern mutex _log_mutex [[maybe_unused]];
#define COUT_DEBUG_FORCE(msg) { scoped_lock<mutex> lock(::_log_mutex); std::cout << "[Stinger::" << __FUNCTION__ << "] [" << common::concurrency::get_thread_id() << "] " << msg << std::endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif


/******************************************************************************
 *                                                                            *
 *  Utility functions                                                         *
 *                                                                            *
 ******************************************************************************/
// dump the content to the given file
static void save(cuckoohash_map<uint64_t, int64_t>& result, const char* dump2file){
    if(dump2file == nullptr) return; // nop
    COUT_DEBUG("save the results to: " << dump2file)

    fstream handle(dump2file, ios_base::out);
    if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

    auto list_entries = result.lock_table();

    for(const auto& p : list_entries){
        handle << p.first << " " << p.second << "\n";
    }

    list_entries.unlock();
    handle.close();
}

/******************************************************************************
 *                                                                            *
 *  Community detection through label propagation                             *
 *                                                                            *
 ******************************************************************************/
namespace gfe::library {

void Stinger::cdlp(uint64_t max_iterations, const char* dump2file){
    TIMER_INIT

    bool change = true;
    int64_t num_mappings = get_max_num_mappings();
    auto ptr_labels0 = make_unique<int64_t[]>(num_mappings); int64_t* labels0 = ptr_labels0.get();
    auto ptr_labels1 = make_unique<int64_t[]>(num_mappings); int64_t* labels1 = ptr_labels1.get();

    // initialisation
    #pragma omp parallel for
    for(int64_t n = 0; n < num_mappings; n++){
        labels0[n] = get_external_id(n);
    }

    // algorithm pass
    for(uint64_t i = 0; i < max_iterations && change; i++){
        COUT_DEBUG("iteration: " << i);
        CHECK_TIMEOUT
        change = false;

        #pragma omp parallel for shared(change)
        for(int64_t n = 0; n < num_mappings; n++){
            labels1[n] = cdlp_propagate(n, labels0);
            if(get_external_id(n) >= 0){
                COUT_DEBUG("label[" << get_external_id(n) << "]  "  << labels0[n] << " -> " << labels1[n] );
            }
            change |= (labels0[n] != labels1[n]);
        }

        std::swap(labels0, labels1); // next iteration
    }

    // save the computation into `dump2file'
    cuckoohash_map<uint64_t, int64_t> result;
    to_external_ids(labels0, num_mappings, &result);
    save(result, dump2file);
}

int64_t Stinger::cdlp_propagate(int64_t vertex_id, int64_t* __restrict labels){
    unordered_map<int64_t, uint64_t> histogram;

    // compute the histogram
    STINGER_FORALL_OUT_EDGES_OF_VTX_BEGIN(STINGER, vertex_id) {
        histogram[labels[STINGER_EDGE_DEST]]++;
    } STINGER_FORALL_OUT_EDGES_OF_VTX_END();

    // cfr. Spec v0.9 pp 14 "If the graph is directed and a neighbor is reachable via both an incoming and
    // outgoing edge, its label will be counted twice"
    if(m_directed){
        STINGER_FORALL_IN_EDGES_OF_VTX_BEGIN(STINGER, vertex_id) {
            histogram[labels[STINGER_EDGE_DEST]]++;
        } STINGER_FORALL_IN_EDGES_OF_VTX_END();
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

    return label_max;
}

} // namespace
