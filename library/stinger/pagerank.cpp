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
#include <fstream>
#include <limits>
#include <mutex>
#include <queue>
#include <thread>

#include "common/system.hpp"
#include "stinger_core/stinger.h"
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
static void save(cuckoohash_map<uint64_t, double>& result, const char* dump2file){
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
 *  Pagerank                                                                  *
 *                                                                            *
 ******************************************************************************/
namespace gfe::library {

// Implementation based on stinger_alg/src/page_rank.c
void Stinger::pagerank(uint64_t num_iterations, double damping_factor, const char* dump2file){
    TIMER_INIT

    const size_t num_registered_vertices = stinger_max_active_vertex(STINGER) +1; // logical vertex
    const size_t num_active_vertices = num_vertices();
    auto ptr_rank0 = make_unique<double[]>(num_registered_vertices); double* rank0t = ptr_rank0.get();
    auto ptr_rank1 = make_unique<double[]>(num_registered_vertices); double* rank1t = ptr_rank1.get();

    // init
    #pragma omp parallel for
    for(uint64_t v = 0; v < num_registered_vertices; v++){
        if(stinger_vtype_get(STINGER, v) == 1) continue; // vertex marked for deletion
        rank0t[v] = 1.0 / num_active_vertices;
    }

    for(uint64_t i = 0; i < num_iterations; i++){
        CHECK_TIMEOUT
        COUT_DEBUG("iteration: " << i);
        double* __restrict rank0 = rank0t; // previous iteration
        double* __restrict rank1 = rank1t; // current iteration

        double pr_constant = 0.0; // dangling sum
        #pragma omp parallel for reduction(+:pr_constant)
        for(uint64_t v = 0; v < num_registered_vertices; v++){
            if(stinger_vtype_get(STINGER, v) == 1) continue; // vertex marked for deletion
            rank1[v] = 0.0;

            // is this vertex a sink? -> add up to the dangling sum
            if(stinger_outdegree(STINGER, v) == 0){ pr_constant += rank0[v]; }

            // compute the score for the current iteration looking at the incoming edges
            STINGER_FORALL_IN_EDGES_OF_VTX_BEGIN(STINGER, v) {
                uint64_t outdegree = stinger_outdegree (STINGER, STINGER_EDGE_DEST);
                if(outdegree > 0){ rank1[v] += rank0[STINGER_EDGE_DEST] / outdegree; }
            } STINGER_FORALL_IN_EDGES_OF_VTX_END();
        }

        #pragma omp parallel for
        for(uint64_t v = 0; v < num_registered_vertices; v++){
            if(stinger_vtype_get(STINGER, v) == 1) continue; // vertex marked for deletion
            rank1[v] = (1.0 - damping_factor) / num_active_vertices + (rank1[v] + (pr_constant / num_active_vertices)) * damping_factor;
        }

        std::swap(rank0t, rank1t); // next iteration
    }


    // store the final results (if required)
    cuckoohash_map<uint64_t, double> result;
    to_external_ids(rank0t, num_registered_vertices, &result); // convert the internal logical IDs into the external IDs
    save(result, dump2file);
};

} // namespace
