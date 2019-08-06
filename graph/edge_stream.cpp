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

#include "edge_stream.hpp"

#include <algorithm>
#include <future>
#include <memory>
#include <thread>
#include <unordered_map>
#include "common/permutation.hpp"
#include "common/timer.hpp"
#include "reader/reader.hpp"
#include "cbytearray.hpp"
#include "configuration.hpp"
#include "edge.hpp"
#include "vertex_list.hpp"

using namespace common;
using namespace std;

namespace graph {

WeightedEdgeStream::WeightedEdgeStream(const std::string& path){
    m_sources = new CByteArray(/* bytes per element */ 8, /* capacity */ 8);
    m_destinations = new CByteArray(/* bytes per element */ 8, /* capacity */ 8);

    auto reader = reader::Reader::open(path);
    WeightedEdge edge;

    Timer timer;
    timer.start();

    while(reader->read(edge)){

        // double the capacity of the arrays sources/destinations/weights
        if(m_num_edges >= m_sources->capacity()){
            auto old_sources = m_sources;
            auto old_destinations = m_destinations;
            auto new_sources = make_unique<CByteArray>(/* bytes per element */ 8, old_sources->capacity() *2);
            auto new_destinations = make_unique<CByteArray>(/* bytes per element */ 8, old_sources->capacity() *2);
            for(size_t i = 0; i < old_sources->capacity(); i++){
                new_sources->set_value_at(i, old_sources->get_value_at(i));
                new_destinations->set_value_at(i, old_destinations->get_value_at(i));
            }

            m_sources = new_sources.release();
            m_destinations = new_destinations.release();

            delete old_sources; old_sources = nullptr;
            delete old_destinations; old_destinations = nullptr;
        }

        m_sources->set_value_at(m_num_edges, edge.m_source);
        m_destinations->set_value_at(m_num_edges, edge.m_destination);
        m_weights.emplace_back(edge.m_weight);

        // keep track of the max vertex id and max weight
        m_max_vertex_id = std::max(m_max_vertex_id, std::max(edge.m_source, edge.m_destination));
        m_max_weight = std::max(m_max_weight, edge.m_weight);

        // update the number of vertices loaded
        m_num_edges++;
    }

    timer.stop();

    LOG("Loaded " << m_num_edges << " edges, max vertex id: " << m_max_vertex_id << ". Load performed in " << timer);
}

WeightedEdgeStream::WeightedEdgeStream(const std::vector<WeightedEdge>& vector){
    m_num_edges = vector.size();
    m_sources = new CByteArray(/* bytes per element */ 8, /* capacity */ m_num_edges);
    m_destinations = new CByteArray(/* bytes per element */ 8, /* capacity */ m_num_edges);
    m_weights.reserve(m_num_edges);

    for(size_t i = 0, sz = m_num_edges; i < sz; i++){
        const auto& edge = vector[i];
        m_sources->set_value_at(i, edge.m_source);
        m_destinations->set_value_at(i, edge.m_destination);
        m_weights.push_back(edge.m_weight);

        m_max_vertex_id = std::max(m_max_vertex_id, std::max(edge.m_source, edge.m_destination));
        m_max_weight = std::max(m_max_weight, edge.m_weight);
    }
}


WeightedEdgeStream::~WeightedEdgeStream(){
    delete m_sources; m_sources = nullptr;
    delete m_destinations; m_destinations = nullptr;
}

void WeightedEdgeStream::permute(){
    permute(configuration().seed() + 91);
}

void WeightedEdgeStream::permute(uint64_t seed){
    if(m_num_edges <= 0) return; // there is nothing to permute

    LOG("Permuting the edge list, seed: " << seed << " ... ");

    Timer timer;
    timer.start();

    // Create a permutation array
    auto ptr_permutation = make_unique<uint64_t[]>(m_num_edges);
    uint64_t* __restrict permutation = ptr_permutation.get();
    for(size_t i = 0; i < m_num_edges; i++){ permutation[i] = i; }
    common::permute(permutation, m_num_edges, seed);

    // Permute the elements in m_sources, m_destinations and m_weights
    auto bytes_per_vertex_id = CByteArray::compute_bytes_per_elements(m_max_vertex_id);
    auto new_sources = make_unique<CByteArray>(/* bytes per element */ bytes_per_vertex_id, m_num_edges);
    auto new_destinations = make_unique<CByteArray>(/* bytes per element */ bytes_per_vertex_id, m_num_edges);
    vector<double> new_weights; new_weights.resize(m_num_edges, 0.0);
//    auto new_weights = make_unique<CByteArray>(/* bytes per element */ CByteArray::compute_bytes_per_elements(m_max_weight), m_num_edges);

    auto permute = [&](uint64_t start, uint64_t length){
        for(size_t i = start, end = start + length; i < end; i++){
            new_sources->set_value_at(i, m_sources->get_value_at(permutation[i]));
            new_destinations->set_value_at(i, m_destinations->get_value_at(permutation[i]));
            new_weights[i] = m_weights[ permutation[i] ];
        }
    };

    uint64_t num_tasks = std::min(m_num_edges, std::max<uint64_t>(4u, thread::hardware_concurrency() * 8));
    uint64_t items_per_task = m_num_edges / num_tasks;
    uint64_t odd_tasks = m_num_edges % num_tasks;
    uint64_t start = 0;
    std::vector<future<void>> tasks;
    tasks.reserve(num_tasks);
    for(size_t i = 0; i < num_tasks; i++){
        uint64_t length = items_per_task + (i < odd_tasks);
        tasks.push_back( async(launch::async, permute, start, length) );
        start += length; // next task
    }
    for(auto& t: tasks) t.get();  // wait for all tasks to finish

    delete m_sources; m_sources = new_sources.release();
    delete m_destinations; m_destinations = new_destinations.release();
    m_weights = std::move(new_weights);

    timer.stop();

    LOG("Permutation completed in " << timer);
}

WeightedEdge WeightedEdgeStream::get(uint64_t index) const {
    if(index >= num_edges()){ INVALID_ARGUMENT("Index out of bound: " << index << " >= " << num_edges()); }
    return WeightedEdge { m_sources->get_value_at(index), m_destinations->get_value_at(index), m_weights[index] };
}

unique_ptr<VertexList> WeightedEdgeStream::vertex_list() const {
    Timer timer;
    timer.start();

    LOG("Computing the list of vertices ...");

    unordered_map<uint64_t, bool> unique_vertices;
    for(uint64_t i = 0, end = m_sources->capacity(); i < end; i++){
        unique_vertices[m_sources->get_value_at(i)] = true;
        unique_vertices[m_destinations->get_value_at(i)] = true;
    }

    auto vertices = make_unique<CByteArray>(CByteArray::compute_bytes_per_elements(m_max_vertex_id), unique_vertices.size());
    uint64_t i = 0;
    for(auto p : unique_vertices){
        vertices->set_value_at(i, p.first);
        i++;
    }

    timer.stop();
    LOG("List of vertices computed in " << timer);

    return make_unique<VertexList>(vertices.release());
}

unique_ptr<cuckoohash_map<uint64_t, bool>> WeightedEdgeStream::vertex_table() const {
    LOG("Computing the list of vertices ... ")
    Timer timer; timer.start();

//    unique_ptr<cuckoohash_map<uint64_t, bool>> ptr_vertex_table;
    auto ptr_vertex_table = make_unique<cuckoohash_map<uint64_t, bool>>();
    auto vertex_table = ptr_vertex_table.get();

    auto populate_vertex_table = [this, vertex_table](uint64_t start, uint64_t length){
        for(uint64_t i = start, end = start + length; i < end; i++){
//            LOG("source: " << m_sources->get_value_at(i) << ", destination: " << m_destinations->get_value_at(i));
            vertex_table->insert(m_sources->get_value_at(i), true);
            vertex_table->insert(m_destinations->get_value_at(i), true);
        }
    };

    const uint64_t num_tasks = thread::hardware_concurrency();
    const uint64_t items_per_task = num_edges() / num_tasks;
    const uint64_t odd_tasks = num_edges() % num_tasks;

    std::vector<future<void>> tasks; tasks.reserve(num_tasks);

    uint64_t start = 0;
    for(size_t i = 0; i < tasks.capacity(); i++){
        uint64_t length = items_per_task + (i < odd_tasks);
        tasks.push_back( async(launch::async, populate_vertex_table, start, length) );
        start += length; // next task
    }
    for(auto& t: tasks) t.get();  // wait for all tasks to finish

    timer.stop();
    LOG("Vertex list computed in " << timer);

    return ptr_vertex_table;
}

}


