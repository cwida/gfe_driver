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
// Stinger stores the weights as int64_t, these macros handle the conversion to double
#define INT2DBL(v) ( *(reinterpret_cast<double*>(&(v))) )
#define DBL2INT(v) ( *(reinterpret_cast<int64_t*>(&(v))) )
#define TIMER_INIT auto time_start = chrono::steady_clock::now();
#define CHECK_TIMEOUT if(timeout.count() > 0 && (chrono::steady_clock::now() - time_start) > timeout) { \
        RAISE_EXCEPTION(TimeoutError, "Timeout occurred after: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - time_start).count() << " seconds") };
#undef CURRENT_ERROR_TYPE
#define CURRENT_ERROR_TYPE ::library::StingerError


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
 *  Library impl~                                                             *
 *                                                                            *
 *****************************************************************************/
namespace library {

// copy & paste from stinger_alg/src/shortest_paths.cpp
namespace {
typedef struct {
    int64_t vertex;
    double  cost;
} weighted_vertex_t;

//this is a simple comparison function for the queue that sorts based on the cost to reach a vertex
bool
static comp(weighted_vertex_t a, weighted_vertex_t b)
{
    return a.cost > b.cost;
}

typedef bool(*CompareType)(weighted_vertex_t a, weighted_vertex_t b);

static std::vector<double> a_star(stinger_t * S, int64_t NV, int64_t source_vertex, int64_t dest_vertex, bool ignore_weights, chrono::seconds timeout){
    TIMER_INIT
    std::vector<double> cost_so_far(NV);

    for (int64_t v = 0; v < NV; v++){
        cost_so_far[v] = std::numeric_limits<double>::max(); // initialize all the values to + infinity
    }

    //represents the vertices to look at in the graph
    std::priority_queue<weighted_vertex_t, std::vector<weighted_vertex_t>, CompareType> frontier (comp);
    //add the initial vertex to the priority queue
    weighted_vertex_t source;
    source.vertex = source_vertex;
    source.cost = 0;
    cost_so_far[source_vertex] = 0;
    frontier.push(source);

    while (!frontier.empty()) {
        CHECK_TIMEOUT

        weighted_vertex_t current = frontier.top();
        frontier.pop();

        //we found our goal!
        if (current.vertex == dest_vertex) break;

        STINGER_FORALL_OUT_EDGES_OF_VTX_BEGIN(S, current.vertex)
                                {
                                    //for all the neighbors of the current vertex
                                    double new_cost;
                                    if (ignore_weights){
                                        new_cost = cost_so_far[current.vertex] + 1; // disregard the weights, instead count paths
                                    } else{
                                        new_cost = cost_so_far[current.vertex] + INT2DBL(STINGER_EDGE_WEIGHT);
                                    }

                                    if (new_cost < cost_so_far[STINGER_EDGE_DEST] || cost_so_far[STINGER_EDGE_DEST] == std::numeric_limits<double>::max()) {
                                        cost_so_far[STINGER_EDGE_DEST] = new_cost;
                                        weighted_vertex_t next;
                                        next.vertex = STINGER_EDGE_DEST;
                                        next.cost = new_cost;
                                        frontier.push(next);
                                    }
                                }
        STINGER_FORALL_OUT_EDGES_OF_VTX_END();
    }

    return cost_so_far;
}



// dump the content to the given file
static void save(cuckoohash_map<uint64_t, double>& result, bool weighted, const char* dump2file){
    if(dump2file == nullptr) return; // nop
    COUT_DEBUG("save the results to: " << dump2file)
    string label_infinity = weighted ? string("infinity") : to_string(numeric_limits<uint64_t>::max());


    fstream handle(dump2file, ios_base::out);
    if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

    auto list_entries = result.lock_table();

    for(const auto& p : list_entries){
        handle << p.first << " ";
        if(p.second == numeric_limits<double>::max()){
            handle << label_infinity;
        } else {
            handle << p.second;
        }
        handle << "\n";
    }

    list_entries.unlock();
    handle.close();
}

} // anonymous namespace


/******************************************************************************
 *                                                                            *
 *  Glue code                                                                 *
 *                                                                            *
 *****************************************************************************/

void Stinger::compute_shortest_paths(uint64_t source_external_id, bool weighted, const char* dump2file){
    int64_t source_internal_id = get_internal_id(source_external_id);
    if(source_internal_id < 0) ERROR("The vertex `" << source_external_id << "' does not exist");
    auto result = a_star(STINGER, stinger_max_active_vertex(STINGER) +1, source_internal_id, /* all vertices */ -1, /* ignore weights */ !weighted, chrono::seconds( m_timeout ));

    // covert the internal vertex IDs back to external IDs
    cuckoohash_map<uint64_t, double> external_ids;
    to_external_ids(result, &external_ids);
    save(external_ids, weighted, dump2file);
}

void Stinger::bfs(uint64_t source_vertex_id, const char* dump2file){
    compute_shortest_paths(source_vertex_id, false, dump2file);
}

void Stinger::sssp(uint64_t source_vertex_id, const char* dump2file){
    compute_shortest_paths(source_vertex_id, true, dump2file);
}

} // namespace
