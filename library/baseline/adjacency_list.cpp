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
#include <iostream>
#include <limits>
#include "reader/reader.hpp"

using namespace std;

namespace library {

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG
static mutex _cout_debug_mutex;
#define COUT_DEBUG_FORCE(msg) { std::scoped_lock<std::mutex> lock{_cout_debug_mutex}; std::cout << "[AdjacencyList::" << __FUNCTION__ << "] " << msg << std::endl; }
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
#define CURRENT_ERROR_TYPE ::library::AdjacencyListError

/*****************************************************************************
 *                                                                           *
 *  Initialisation                                                           *
 *                                                                           *
 *****************************************************************************/

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
    scoped_lock<mutex_t> lock(m_mutex);
    return m_adjacency_list.size();
}

uint64_t AdjacencyList::num_edges() const {
    scoped_lock<mutex_t> lock(m_mutex);
    return m_num_edges;
}

bool AdjacencyList::has_vertex(uint64_t vertex_id) const {
    scoped_lock<mutex_t> lock(m_mutex);
    return m_adjacency_list.find(vertex_id) != end(m_adjacency_list);
}

double AdjacencyList::get_weight(uint64_t source, uint64_t destination) const {
    constexpr double NaN { numeric_limits<double>::signaling_NaN() };

    scoped_lock<mutex_t> lock(m_mutex);
    auto vertex_src = m_adjacency_list.find(source);
    if(vertex_src == end(m_adjacency_list)) return NaN;

    auto& outgoing_edges = vertex_src->second.first;
    auto result = find_if(begin(outgoing_edges), end(outgoing_edges), [destination](const pair<uint64_t, double>& edge){
        return edge.first == destination;
    });

    if(result == end(outgoing_edges)) return NaN;
    return result->second;
}

/*****************************************************************************
 *                                                                           *
 *  Updates                                                                  *
 *                                                                           *
 *****************************************************************************/

bool AdjacencyList::add_vertex(uint64_t vertex_id){
    scoped_lock<mutex_t> lock(m_mutex);
    COUT_DEBUG("vertex_id: " << vertex_id);
    auto pair = m_adjacency_list.emplace( vertex_id, EdgePair{} );
    return pair.second;
}

bool AdjacencyList::delete_vertex(uint64_t vertex_id){
    scoped_lock<mutex_t> lock(m_mutex);
    COUT_DEBUG("vertex_id: " << vertex_id);
    auto vertex_src = m_adjacency_list.find(vertex_id);
    if(vertex_src == end(m_adjacency_list)) return false;

    // remove all dangling edges
    for(auto& e : vertex_src->second.first){
        COUT_DEBUG("delete edge " << vertex_src << " -> " << e.first);
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
            COUT_DEBUG("delete edge " << vertex_src << " <- " << e.first);
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
    scoped_lock<mutex_t> lock(m_mutex);
    COUT_DEBUG("edge: " << e);

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

bool AdjacencyList::delete_edge(graph::Edge e){
    scoped_lock<mutex_t> lock(m_mutex);
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
        add_vertex(edge.m_source);
        add_vertex(edge.m_destination);
        add_edge(edge);
    }
}


/*****************************************************************************
 *                                                                           *
 *  Dump                                                                     *
 *                                                                           *
 *****************************************************************************/
void AdjacencyList::dump(std::ostream& out) const {
    scoped_lock<mutex_t> lock(m_mutex);

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
