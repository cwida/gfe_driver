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
#include <cinttypes>
#include <limits>
#include <mutex>
#include <queue>
#include <thread>

#include "common/system.hpp"
#include "stinger_core/stinger.h"

using namespace std;

namespace library {

// Macros
#define STINGER reinterpret_cast<struct stinger*>(m_stinger_graph)
#undef CURRENT_ERROR_TYPE
#define CURRENT_ERROR_TYPE ::library::StingerError

/******************************************************************************
 *                                                                            *
 *  Debug                                                                     *
 *                                                                            *
 *****************************************************************************/
extern mutex _stinger_debug_mutex;
//#define DEBUG
#define COUT_DEBUG_FORCE(msg) { scoped_lock<mutex> lock(_stinger_debug_mutex); std::cout << "[Stinger::" << __FUNCTION__ << "] [" << common::concurrency::get_thread_id() << "] " << msg << std::endl; }
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
// copy & paste from stinger_alg/src/shortest_paths.cpp
namespace {
typedef struct{
    int64_t vertex;
    int64_t cost;
}weighted_vertex_t;

//this is a simple comparison function for the queue that sorts based on the cost to reach a vertex
bool
comp(weighted_vertex_t a, weighted_vertex_t b)
{
    return a.cost > b.cost;
}
typedef bool(*CompareType)(weighted_vertex_t a, weighted_vertex_t b);

static int64_t a_star(stinger_t * S, int64_t NV, int64_t source_vertex, int64_t dest_vertex, bool ignore_weights){
    std::vector<int64_t> cost_so_far(NV);
    for (int64_t v = 0; v < NV; v++){
        cost_so_far[v] = std::numeric_limits<int64_t>::max(); // initialize all the values to infinity
    }
    //represents the verties to look at in the graph
    std::priority_queue<weighted_vertex_t, std::vector<weighted_vertex_t>, CompareType> frontier (comp);
    //add the initial vertex to the priority queue
    weighted_vertex_t source;
    source.vertex = source_vertex;
    source.cost = 0;
    cost_so_far[source_vertex] = 0;
    frontier.push(source);

    while (!frontier.empty()) {
        weighted_vertex_t current = frontier.top();
        frontier.pop();

        if (current.vertex == dest_vertex) {
            //we found our goal!
            return cost_so_far[current.vertex];
        }

        STINGER_FORALL_OUT_EDGES_OF_VTX_BEGIN(S, current.vertex)
                                {
                                    //for all the neighbors of the current vertex
                                    int64_t new_cost;
                                    if (ignore_weights){
                                        new_cost= cost_so_far[current.vertex] + 1; // disregard the weights, instead count paths
                                    }
                                    else{
                                        new_cost = cost_so_far[current.vertex] + STINGER_EDGE_WEIGHT;
                                    }

                                    if (new_cost < cost_so_far[STINGER_EDGE_DEST] || cost_so_far[STINGER_EDGE_DEST] == std::numeric_limits<int64_t>::max()) {
                                        cost_so_far[STINGER_EDGE_DEST] = new_cost;
                                        weighted_vertex_t next;
                                        next.vertex = STINGER_EDGE_DEST;
                                        next.cost = new_cost;
                                        frontier.push(next);
                                    }
                                }
        STINGER_FORALL_OUT_EDGES_OF_VTX_END();

    }

    //we failed to reach this vertex - meaning that there is no path between these two verticies
    if(dest_vertex >= 0)
        return cost_so_far[dest_vertex]; //should return max integer value as a placeholder for infinity
    else
        return 0; // shortest path towards all vertices
}

static int64_t a_star_debug(stinger_t * S, int64_t NV, int64_t source_vertex, int64_t dest_vertex, bool ignore_weights, std::vector<library::ShortestPathInterface::Distance>& result){
    std::vector<int64_t> cost_so_far(NV);
    std::vector<int64_t> predecessor(NV);
    for (int64_t v = 0; v < NV; v++){
        cost_so_far[v] = std::numeric_limits<int64_t>::max(); // initialize all the values to infinity
        predecessor[v] = -1;
    }
    //represents the verties to look at in the graph
    std::priority_queue<weighted_vertex_t, std::vector<weighted_vertex_t>, CompareType> frontier (comp);
    //add the initial vertex to the priority queue
    weighted_vertex_t source;
    source.vertex = source_vertex;
    source.cost = 0;
    cost_so_far[source_vertex] = 0;
    predecessor[source_vertex] = source_vertex;
    frontier.push(source);

    while (!frontier.empty()) {
        weighted_vertex_t current = frontier.top();
        frontier.pop();

        if (current.vertex == dest_vertex) {
            break; //we found our goal!
        } else {
            STINGER_FORALL_OUT_EDGES_OF_VTX_BEGIN(S, current.vertex)
                                    {
                                        //for all the neighbors of the current vertex
                                        int64_t new_cost;
                                        if (ignore_weights){
                                            new_cost= cost_so_far[current.vertex] + 1; // disregard the weights, instead count paths
                                        }
                                        else{
                                            new_cost = cost_so_far[current.vertex] + STINGER_EDGE_WEIGHT;
                                        }

                                        if (new_cost < cost_so_far[STINGER_EDGE_DEST] || cost_so_far[STINGER_EDGE_DEST] == std::numeric_limits<int64_t>::max()) {
                                            cost_so_far[STINGER_EDGE_DEST] = new_cost;
                                            predecessor[STINGER_EDGE_DEST] = STINGER_EDGE_SOURCE;

                                            weighted_vertex_t next;
                                            next.vertex = STINGER_EDGE_DEST;
                                            next.cost = new_cost;
                                            frontier.push(next);
                                        }
                                    }
            STINGER_FORALL_OUT_EDGES_OF_VTX_END();
        }
    }

    result.clear();

    if(dest_vertex >= 0) {
        if(cost_so_far[dest_vertex] < std::numeric_limits<int64_t>::max()){
            int64_t v = predecessor[dest_vertex];
            while(v != source_vertex){
                result.emplace_back(v, cost_so_far[v]);
                v = predecessor[v];
            }
            std::reverse(begin(result), end(result));
        }

        return cost_so_far[dest_vertex]; //should return max integer value as a placeholder for infinity
    } else {
        for (int64_t v = 0; v < NV; v++){
            if(cost_so_far[v] < std::numeric_limits<int64_t>::max()){
                result.emplace_back(v, cost_so_far[v]);
            }
        }
        return 0; // shortest path towards all vertices
    }

}


} // anonymouse namespace


/******************************************************************************
 *                                                                            *
 *  Glue code                                                                 *
 *                                                                            *
 *****************************************************************************/

int64_t Stinger::compute_shortest_paths(uint64_t ext_source, uint64_t ext_destination, bool weighted, std::vector<library::ShortestPathInterface::Distance>* result){
    // internal ids in the adjacency list
    int64_t id_source = get_internal_id(ext_source);
    if(id_source < 0){
        COUT_DEBUG("The source vertex does not exist: " << ext_source);
        return -1;
    }
    int64_t id_destination = -1;
    if(id_destination != std::numeric_limits<uint64_t>::max()){
        id_destination = get_internal_id(ext_destination);
        if(id_destination < 0){
            COUT_DEBUG("The destination vertex does not exist: " << ext_source);
            return -1;
        }
    }

    // trampoline
    int64_t distance {0};
    if(result == nullptr){
        distance = a_star(STINGER, stinger_max_active_vertex(STINGER), id_source, id_destination, !weighted);
    } else {
        distance = a_star_debug(STINGER, stinger_max_active_vertex(STINGER), id_source, id_destination, !weighted, *result);

        // transform the internal IDs back to the external vertex IDs
        for(uint64_t i =0, sz = result->size(); i < sz; i++){
            (*result)[i].m_vertex = get_external_id( (*result)[i].m_vertex );
        }
    }

    return distance;
}

void Stinger::bfs_all(uint64_t source, std::vector<library::ShortestPathInterface::Distance>* result){
    compute_shortest_paths(source, std::numeric_limits<uint64_t>::max(), false, result);
}

int64_t Stinger::bfs_one(uint64_t source, uint64_t dest, std::vector<library::ShortestPathInterface::Distance>* path){
    return compute_shortest_paths(source, dest, false, path);
}

void Stinger::spw_all(uint64_t source, std::vector<library::ShortestPathInterface::Distance>* result){
    compute_shortest_paths(source, std::numeric_limits<uint64_t>::max(), true, result);
}

int64_t Stinger::spw_one(uint64_t source, uint64_t dest, std::vector<library::ShortestPathInterface::Distance>* path){
    return compute_shortest_paths(source, dest, true, path);
}

} // namespace
