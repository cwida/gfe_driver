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

#include "adjacency_list.hpp"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <mutex>
#include <queue>
#include <unordered_set>

#include "common/circular_array.hpp"
#include "common/system.hpp"
#include "common/timer.hpp"
#include "configuration.hpp" // LOG
#include "reader/reader.hpp"
#include "third-party/robin_hood/robin_hood.h"

using namespace common;
using namespace std;

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG
extern mutex _log_mutex;
#define COUT_DEBUG_FORCE(msg) { std::scoped_lock<std::mutex> lock{_log_mutex}; std::cout << "[AdjacencyList::" << __FUNCTION__ << "] [" << concurrency::get_thread_id() << "] " << msg << std::endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif

/*****************************************************************************
 *                                                                           *
 *  Error                                                                    *
 *                                                                           *
 *****************************************************************************/
#undef CURRENT_ERROR_TYPE
#define CURRENT_ERROR_TYPE ::gfe::library::AdjacencyListError

/*****************************************************************************
 *                                                                           *
 *  Initialisation                                                           *
 *                                                                           *
 *****************************************************************************/
namespace gfe::library {


AdjacencyList::AdjacencyList(bool is_directed) : m_is_directed(is_directed) {

}

AdjacencyList::~AdjacencyList(){ }

/*****************************************************************************
 *                                                                           *
 *  Properties                                                               *
 *                                                                           *
 *****************************************************************************/

bool AdjacencyList::is_directed() const {
    // m_is_directed is const, no need for a critical section here
    return m_is_directed;
}

uint64_t AdjacencyList::num_vertices() const {
    shared_lock<mutex_t> lock(m_mutex);
    return m_adjacency_list.size();
}

uint64_t AdjacencyList::num_edges() const {
    shared_lock<mutex_t> lock(m_mutex);
    return m_num_edges;
}

bool AdjacencyList::has_vertex(uint64_t vertex_id) const {
    shared_lock<mutex_t> lock(m_mutex);
    return m_adjacency_list.find(vertex_id) != end(m_adjacency_list);
}

double AdjacencyList::get_weight(uint64_t source, uint64_t destination) const {
    constexpr double NaN { numeric_limits<double>::signaling_NaN() };

    shared_lock<mutex_t> lock(m_mutex);
    auto vertex_src = m_adjacency_list.find(source);
    if(vertex_src == end(m_adjacency_list)) return NaN;

    auto& outgoing_edges = vertex_src->second.first;
    auto result = find_if(begin(outgoing_edges), end(outgoing_edges), [destination](const pair<uint64_t, double>& edge){
        return edge.first == destination;
    });

    if(result == end(outgoing_edges)) return NaN;
    return result->second;
}

uint64_t AdjacencyList::get_degree(uint64_t vertex_id) const {
    auto vertex = m_adjacency_list.find(vertex_id);
    assert(vertex != end(m_adjacency_list) && "The given edge does not exist");
    if(!is_directed()){ // undirected
        return vertex->second.first.size();
    } else {
        std::unordered_set<uint64_t> neighbours;
        neighbours.reserve(vertex->second.first.size() + vertex->second.second.size());
        for(const auto& e: vertex->second.first){ neighbours.insert(e.first); }
        for(const auto& e: vertex->second.second){ neighbours.insert(e.first); }
        return neighbours.size();
    }
}

const AdjacencyList::EdgeList& AdjacencyList::get_incoming_edges(uint64_t vertex_id) const {
    auto it = m_adjacency_list.find(vertex_id);
    if(it == end(m_adjacency_list)) ERROR("The searched vertex `" << vertex_id << "' does not exist");
    return get_incoming_edges(it->second);
}

const AdjacencyList::EdgeList& AdjacencyList::get_incoming_edges(const NodeList::mapped_type& adjlist_value) const {
    return is_directed() ? /* directed graph */ adjlist_value.second : /* undirected graph */ adjlist_value.first;
}

bool AdjacencyList::has_timeout() const {
    return m_timeout.count() > 0;
}

void AdjacencyList::set_timeout(uint64_t seconds){
    COUT_DEBUG("Timeout set to " << seconds << " seconds");
    m_timeout = chrono::seconds{seconds};
}

/*****************************************************************************
 *                                                                           *
 *  Helpers                                                                  *
 *                                                                           *
 *****************************************************************************/
template<typename Function>
void AdjacencyList::for_all_edges(const EdgePair& p, Function f){
    if(!is_directed()){
        for(const auto& e : p.first){ std::invoke(f, e.first); }
    } else {
        // avoid invoking f on the same edge having the same two endpoints ( a -> b and a <- b )
        std::unordered_set<uint64_t> neighbours;
        neighbours.reserve(p.first.size() + p.second.size());
        for(const auto& e : p.first){ neighbours.insert(e.first); }
        for(const auto& e : p.second){ neighbours.insert(e.first); }
        for(const auto& v : neighbours) { std::invoke(f, v); }
    }
}

bool AdjacencyList::check_directed_edge_exists(uint64_t vertex1, uint64_t vertex2){
    const auto& p = m_adjacency_list.find(vertex1);
    assert(p != end(m_adjacency_list) && "The node `vertex1' does not exist");
    assert(m_adjacency_list.find(vertex2) != end(m_adjacency_list) && "The node `vertex2' does not exist");

    for(const auto& e: p->second.first){
        if(e.first == vertex2) return true;
    }

    return false;
}

/*****************************************************************************
 *                                                                           *
 *  Updates                                                                  *
 *                                                                           *
 *****************************************************************************/

bool AdjacencyList::add_vertex(uint64_t vertex_id){
    scoped_lock<mutex_t> lock(m_mutex);
    return add_vertex0(vertex_id);
}

bool AdjacencyList::add_vertex0(uint64_t vertex_id){
    COUT_DEBUG("vertex_id: " << vertex_id);
    auto pair = m_adjacency_list.emplace( vertex_id, EdgePair{} );
    return pair.second;
}

bool AdjacencyList::remove_vertex(uint64_t vertex_id){
    scoped_lock<mutex_t> lock(m_mutex);
    return delete_vertex0(vertex_id);
}

bool AdjacencyList::delete_vertex0(uint64_t vertex_id){
    COUT_DEBUG("vertex_id: " << vertex_id);
    auto vertex_src = m_adjacency_list.find(vertex_id);
    if(vertex_src == end(m_adjacency_list)) return false;

    // remove all dangling edges
    for(auto& e : vertex_src->second.first){
        COUT_DEBUG("delete edge " << vertex_src->first << " -> " << e.first);
        auto vertex_dst = m_adjacency_list.find(e.first);
        assert(vertex_dst != end(m_adjacency_list) && "the dest vertex does not exist (#1)");
        auto lst_back = m_is_directed ? vertex_dst->second.second : vertex_dst->second.first;
        auto it = find_if(begin(lst_back), end(lst_back), [vertex_id](const pair<uint64_t,double>& edge){
            return edge.first == vertex_id;
        });
        assert(it != end(lst_back) && "dangling edge");
        lst_back.erase(it);
    }
    assert(m_is_directed || vertex_src->second.second.empty()); // if the graph is undirected, then the list of incoming edges should be empty!
    if(m_is_directed){
        for(auto& e : vertex_src->second.second){
            COUT_DEBUG("delete edge " << vertex_src->first << " <- " << e.first);
            auto vertex_dst = m_adjacency_list.find(e.first);
            assert(vertex_dst != end(m_adjacency_list) && "the dest vertex does not exist (#2)");
            auto lst_fwd = vertex_dst->second.first;
            auto it = find_if(begin(lst_fwd), end(lst_fwd), [vertex_id](const pair<uint64_t,double>& edge){
                return edge.first == vertex_id;
            });
            assert(it != end(lst_fwd) && "dangling edge");
            lst_fwd.erase(it);
        }
    }

    m_adjacency_list.erase(vertex_src);

    return true;
}

bool AdjacencyList::add_edge(graph::WeightedEdge e){
    if(e.source() == e.destination()) INVALID_ARGUMENT("Cannot insert an edge with the same source and destination: " << e);

    scoped_lock<mutex_t> lock(m_mutex);
    return add_edge0(e);
}

bool AdjacencyList::add_edge0(graph::WeightedEdge e){
    COUT_DEBUG("edge: " << e);
    assert(e.m_weight >= 0 && "Negative weight");

    auto v_src = m_adjacency_list.find(e.source());
    if(v_src == m_adjacency_list.end()){
        COUT_DEBUG("The source vertex " << e.source() << " does not exist");
        return false;
    }
    auto v_dst = m_adjacency_list.find(e.destination());
    if(v_dst == m_adjacency_list.end()){
        COUT_DEBUG("The destination vertex " << e.destination() << " does not exist");
        return false;
    }

    auto& list_out = v_src->second.first;
    auto it = find_if(begin(list_out), end(list_out), [e](const pair<uint64_t,double>& edge){
       return edge.first == e.destination();
    });
    if(it == end(list_out)){
        list_out.emplace_back(e.destination(), e.weight());
        m_num_edges++;
    } else {
        it->second = e.weight();
    }

    auto& list_in = m_is_directed ? v_dst->second.second : v_dst->second.first;
    it = find_if(begin(list_in), end(list_in), [e](const pair<uint64_t, double>& edge){
       return e.source() == edge.first;
    });
    if(it == end(list_in)){
        list_in.emplace_back(e.source(), e.weight());
    } else {
        it->second = e.weight();
    }

    return true;
}

bool AdjacencyList::remove_edge(graph::Edge e){
    scoped_lock<mutex_t> lock(m_mutex);
    return delete_edge0(e);
}

bool AdjacencyList::delete_edge0(graph::Edge e){
    COUT_DEBUG("edge: " << e);

    auto vertex_src = m_adjacency_list.find(e.source());
    if(vertex_src == end(m_adjacency_list)) return false;
    auto vertex_dst = m_adjacency_list.find(e.destination());
    if(vertex_dst == end(m_adjacency_list)) return false;

    auto& list_out = vertex_src->second.first;
    auto it = find_if(begin(list_out), end(list_out), [e](const pair<uint64_t, double>& edge){
        return e.destination() == edge.first;
    });
    if(it != end(list_out)){
        list_out.erase(it);
        assert(m_num_edges > 0 && "underflow");
        m_num_edges--;
    } else {
        return false; // if there is no outgoing edge from source to destination, there cannot be an incoming edge from destination to source
    }

    auto& list_in = m_is_directed ? vertex_dst->second.second : vertex_dst->second.first;
    it = find_if(begin(list_in), end(list_in), [e](const pair<uint64_t, double>& edge){
        return e.source() == edge.first;
    });
    assert(it != end(list_in) && "the outgoing edge was present, but no incoming edge");
    list_in.erase(it);

    return true;
}

void AdjacencyList::load(const std::string& path){
    scoped_lock<mutex_t> lock(m_mutex);
    COUT_DEBUG("path: " << path);
    auto reader = reader::Reader::open(path);
    ASSERT(reader->is_directed() == m_is_directed);
    graph::WeightedEdge edge;
    while(reader->read(edge)){
        add_vertex0(edge.m_source);
        add_vertex0(edge.m_destination);
        add_edge0(edge);
    }
}


/*****************************************************************************
 *                                                                           *
 *  Graphalytics                                                             *
 *                                                                           *
 *****************************************************************************/
#define TIMER_INIT auto time_start = chrono::steady_clock::now();
#define CHECK_TIMEOUT if(has_timeout() && (chrono::steady_clock::now() - time_start) > m_timeout) { \
        RAISE_EXCEPTION(TimeoutError, "Timeout occurred after: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - time_start).count() << " seconds") };

void AdjacencyList::bfs(uint64_t source_vertex_id, const char* dump2file){
    TIMER_INIT
    shared_lock<mutex_t> lock(mutex);

    // init
    unordered_map<uint64_t, int64_t> distances; distances[source_vertex_id] = 0;
    common::CircularArray<uint64_t> queue; queue.append(source_vertex_id);

    while(!queue.empty()){
        CHECK_TIMEOUT

        uint64_t vertex = queue[0]; queue.pop();
        int64_t distance = distances.at(vertex);
        const auto& outgoing_edges = m_adjacency_list.at(vertex).first;
        for(auto& edge : outgoing_edges){
            auto result = distances.emplace(edge.first, distance +1);
            if(result.second){ queue.append(edge.first); }
        }
    }

    if(dump2file != nullptr){
        COUT_DEBUG("save the results to: " << dump2file)
        fstream handle(dump2file, ios_base::out);
        if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

        for(const auto& v : m_adjacency_list){
            handle << v.first << " ";
            auto distance = distances.find(v.first);
            if(distance != end(distances)){
                handle << distance->second;
            } else {
                handle << std::numeric_limits<int64_t>::max(); // it should have been -1, but ok
            }
            handle << "\n";
        }

        handle.close();
    }
}

void AdjacencyList::pagerank(uint64_t num_iterations, double damping_factor, const char* dump2file){
    TIMER_INIT
    shared_lock<mutex_t> lock(mutex);

    // init
    unordered_map<uint64_t, double> rank;
    unordered_map<uint64_t, double> rank_next;
    unordered_set<uint64_t> sinks; // vertices with no outgoing edges
    const uint64_t V = num_vertices();
    for(const auto& p : m_adjacency_list){
        rank[p.first] = 1.0/V;
        if(p.second.first.empty()) sinks.insert(p.first);
    }

    // perform `num_iterations' of the pagerank algorithm
    for(uint64_t iter = 0; iter < num_iterations; iter++){
        CHECK_TIMEOUT

        // step #1, compute the `leakage', the score repartitioned from the sinks to all the other nodes in the graph
        // lines 5 - 10 of the spec v1.0 pp 36
        double dangling_sum = 0;
        for(uint64_t v : sinks){ dangling_sum += rank[v]; }

        // step #2, compute the rank for the current iteration, for all vertices
        for(const auto& p : m_adjacency_list){
            uint64_t vertex_id = p.first;
            double score = 0;
            const auto& incoming_edges = get_incoming_edges(p.second);
            for(const auto& e: incoming_edges){
//                uint64_t source_id = e.first == vertex_id ? e.second : e.first; // in case the graph is undirected, the vertex id may be reversed
                uint64_t source_id = e.first;
                uint64_t nout = m_adjacency_list.find(source_id)->second.first.size();
                if(nout > 0){ score += rank[source_id] / nout; }
            }
            // final score for the given vertex_id
            score = (1.0 - damping_factor) / V + damping_factor * score + (damping_factor/V) * dangling_sum;
            rank_next[vertex_id] = score;
        }

        // set the computed rank the current for the next iteration
        std::swap(rank, rank_next);
    }

    if(dump2file != nullptr){
        COUT_DEBUG("save the results to: " << dump2file)
        fstream handle(dump2file, ios_base::out);
        if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

        for(const auto& v : rank){
            handle << v.first << " " << v.second << "\n";
        }

        handle.close();
    }

}

void AdjacencyList::wcc(const char* dump2file){
    TIMER_INIT
    shared_lock<mutex_t> lock(mutex);

    // init
    unordered_map<uint64_t, uint64_t> components;
    for(const auto& p: m_adjacency_list) components[p.first] = p.first;

    // repeat until it converges
    bool converged {false};
    do {
        CHECK_TIMEOUT
        converged = true;

        for(const auto& p: m_adjacency_list){
            uint64_t vertex_id = p.first;
            uint64_t min_component = numeric_limits<int64_t>::max();

            // outgoing edges
            for(const auto& out : p.second.first){
                min_component = min(min_component, components[out.first]);
            }
            // incoming edges
            for(const auto& in : p.second.second){
                min_component = min(min_component, components[in.first]);
            }

            // update?
            auto it_component = components.find(vertex_id);
            assert(it_component != end(components) && "vertex_id not found");
            if(it_component->second > min_component){
                it_component->second = min_component;
                converged = false; // and repeat
            }
        }
    } while(!converged);

    if(dump2file != nullptr){
        COUT_DEBUG("save the results to: " << dump2file)
        fstream handle(dump2file, ios_base::out);
        if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

        for(const auto& c : components){
            handle << c.first << " " << c.second << "\n";
        }

        handle.close();
    }
}

void AdjacencyList::cdlp(uint64_t max_iterations, const char* dump2file){
    TIMER_INIT
    shared_lock<mutex_t> lock(m_mutex);

    // init
    unordered_map<uint64_t, uint64_t> labels;
    unordered_map<uint64_t, uint64_t> labels_next;
    unordered_map<uint64_t, uint64_t> histogram;

    // first iteration
    for(auto& v : m_adjacency_list){
        labels[v.first] = v.first;
    }

    // bulk of the algorithm, perform up to `max_iterations'
    uint64_t iteration = 1;
    bool change = true;
    while(iteration <= max_iterations && change){
        COUT_DEBUG("iteration: " << iteration);

        CHECK_TIMEOUT
        change = false;

        // for each vertex...
        for(auto& v : m_adjacency_list){
            histogram.clear();

            // count the labels among the neighbours of v
            for(auto& e : v.second.first){ // outgoing edges
                histogram[ labels[e.first] ] += 1;
            }
            for(auto& e : v.second.second){ // incoming edges
                histogram[ labels[e.first] ] += 1;
            }

            // get the label with the max frequency && the smallest id
            uint64_t label_id = 0, max_count = 0;
            for(auto& h: histogram){
                if(h.second > max_count) {
                    max_count = h.second;
                    label_id = h.first;
                } else if (h.second == max_count && h.first < label_id){ // select the minimum label with the highest count
                    label_id = h.first;
                }
            }
//            COUT_DEBUG("label selected: " << label_id << ", count: " << max_count);

            labels_next[v.first] = label_id;
            COUT_DEBUG("label[" << v.first << "]: " << labels[v.first] << " -> " << labels_next[v.first]);
            change |= ( labels[v.first] != labels_next[v.first] ); // did we update the label ?
        }

        // prepare for the next iteration
        swap(labels, labels_next);
        iteration++;
    }

    if(dump2file != nullptr){
        COUT_DEBUG("save the results to: " << dump2file)
        fstream handle(dump2file, ios_base::out);
        if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

        for(auto& l : labels){
            handle << l.first << " " << l.second << "\n";
        }

        handle.close();
    }
}

void AdjacencyList::lcc_undirected(unordered_map<uint64_t, double>& result){
    TIMER_INIT
    robin_hood::unordered_map<uint64_t, bool> H;
    robin_hood::unordered_map<uint64_t, uint64_t> already_visited;

    vector<uint64_t> nodes; nodes.reserve(m_adjacency_list.size());
    for(const auto& p: m_adjacency_list){ nodes.push_back(p.first); }
    std::sort(begin(nodes), end(nodes));

    Timer t_build, t_probe;
    uint64_t hash_probes = 0;

    for(auto vertex1 : nodes){
        CHECK_TIMEOUT
        uint64_t degree = get_degree(vertex1);

        if(degree >= 2){
            H.clear();
            const EdgePair& vertex1_neighbours = m_adjacency_list[vertex1];
            uint64_t num_triangles = already_visited[vertex1];
//            COUT_DEBUG("vertex1: " << vertex1 << ", degree: " << degree << ", num_triangles: " << num_triangles);

            // build the hash table
            t_build.resume();
            for(auto& p: vertex1_neighbours.first){ H[p.first] = true; }
            t_build.stop();

            // probe the hash table
            t_probe.resume();
            for(auto& p: vertex1_neighbours.first){
                uint64_t vertex2 = p.first;
                if(vertex2 < vertex1) continue; // already visited
//                COUT_DEBUG("vertex1: " << vertex1 << ", vertex2: " << vertex2);

                const EdgeList& vertex2_neighbours = m_adjacency_list.at(vertex2).first;
                for(auto& edge : vertex2_neighbours){
                    hash_probes ++;
                    bool triangle_found = H.find(edge.first) != cend(H);
                    if(triangle_found){
                        COUT_DEBUG("triangle found: " << vertex1 << " - " << vertex2 << " - " << edge.first);
                        num_triangles++;
                        already_visited[vertex2] ++;
                    }
                }
            }

            t_probe.stop();

            result[vertex1] = static_cast<double>(num_triangles) / (degree * (degree-1));
        } else {
            result[vertex1] = 0;
        }

    }

    LOG("timer build: " << t_build << ", timer probe: " << t_probe << ", number hashes: " << hash_probes);
}

void AdjacencyList::lcc_directed(unordered_map<uint64_t, double>& result){
    TIMER_INIT
    for(const auto& p : m_adjacency_list){
        CHECK_TIMEOUT
        uint64_t vertex1 = p.first;
        uint64_t degree = get_degree(vertex1);

        if(degree >= 2){
            uint64_t num_triangles = 0;

            for_all_edges(p.second, [&](uint64_t vertex2){
                for_all_edges(p.second, [&, vertex2](uint64_t vertex3){
                    if(vertex2 != vertex3){
#if defined(DEBUG)
                        if( check_directed_edge_exists(vertex2, vertex3) ){
                            COUT_DEBUG("triangle found: " << vertex1 << " - " << vertex2 << " - " << vertex3);
                        }
#endif
                        num_triangles += check_directed_edge_exists(vertex2, vertex3);
                    }
                });
            });

            result[vertex1] = static_cast<double>(num_triangles) / (degree * (degree-1));
        } else {
            result[vertex1] = 0;
        }
    }
}

void AdjacencyList::lcc(const char* dump2file){
    shared_lock<mutex_t> lock(mutex);

    // init
    std::unordered_map<uint64_t, double> lcc;

    // computation
    if(is_directed()){
        lcc_directed(lcc);
    } else {
        lcc_undirected(lcc);
    }

    if(dump2file != nullptr){
        COUT_DEBUG("save the results to: " << dump2file)
        fstream handle(dump2file, ios_base::out);
        if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

        for(const auto& p : lcc){
            handle << p.first << " " << p.second << "\n";
        }

        handle.close();
    }
}

void AdjacencyList::sssp(uint64_t source_vertex_id, const char* dump2file){
    TIMER_INIT
    shared_lock<mutex_t> lock(mutex);

    // init
    std::unordered_set<uint64_t> visited;
    std::unordered_map<uint64_t, double> distances;
    for(auto& v : m_adjacency_list){ distances[v.first] = numeric_limits<double>::infinity(); }
    distances[source_vertex_id] = 0.0;

    // priority queue
    struct elt { uint64_t m_vertex_id; double m_distance; };
    auto comparator = [](const elt& e1, const elt& e2){ return e1.m_distance > e2.m_distance; };
    priority_queue<elt, vector<elt>, decltype(comparator)> Q(comparator);
    Q.push({source_vertex_id, 0.0});

    while(!Q.empty()){ // Dijkstra
        CHECK_TIMEOUT
        auto elt = Q.top(); Q.pop();
        COUT_DEBUG("extract: " << elt.m_vertex_id << ", distance: " << elt.m_distance);
        if( !visited.insert(elt.m_vertex_id).second ) continue; // we already processed this vertex

        const auto& out_edges = m_adjacency_list[elt.m_vertex_id].first;
        for(const auto& e : out_edges){
            if(elt.m_distance + e.second < distances[e.first]){
                distances[e.first] = elt.m_distance + e.second;
                COUT_DEBUG("push: " << e.first << ", distance: " << elt.m_distance + e.second);
                Q.push({ e.first, elt.m_distance + e.second } );
            }
        }
    }

    if(dump2file != nullptr){
        COUT_DEBUG("save the results to: " << dump2file)
        fstream handle(dump2file, ios_base::out);
        if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

        for(const auto& d : distances){
            handle << d.first << " " << d.second << "\n";
        }

        handle.close();
    }
}

/*****************************************************************************
 *                                                                           *
 *  Dump                                                                     *
 *                                                                           *
 *****************************************************************************/
void AdjacencyList::dump_ostream(std::ostream& out) const {
    shared_lock<mutex_t> lock(m_mutex);

    out << "[AdjacencyList] vertices: " << num_vertices() << ", edges: " << num_edges() << ", "
            "graph " << (is_directed() ? "directed" : "undirected") << "\n";

    vector<uint64_t> vertices; vertices.reserve(num_vertices());
    for(auto& p : m_adjacency_list){ vertices.push_back(p.first); }
    sort(begin(vertices), end(vertices));

    for(uint64_t vertex_src_id : vertices){
        out << "[" << vertex_src_id << "] outgoing edges: ";
        const auto& outgoing_edges = m_adjacency_list.at(vertex_src_id).first;
        bool first = true;
        for(const auto& edge : outgoing_edges){
            if(first){ first = false; } else { out << ", "; }
            out << edge.first << " (" << edge.second << ")";
        }
        out << "\n";
        if(is_directed()){
            out << "\tincoming edges: ";
            const auto& incoming_edges = m_adjacency_list.at(vertex_src_id).second;

            first = true;
            for(const auto& edge : incoming_edges){
                if(first){ first = false; } else { out << ", "; }
                out << edge.first << " (" << edge.second << ")";
            }
            out << "\n";
        }
    }
}


} // namespace
