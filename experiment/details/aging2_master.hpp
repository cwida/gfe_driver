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

#pragma once

#include <memory>
#include <vector>

#include "common/static_index.hpp"
#include "experiment/aging2_result.hpp"
#include "graph/edge.hpp"
#include "third-party/libcuckoo/cuckoohash_map.hh"

// forward declarations
namespace experiment { class Aging2Experiment; }
namespace experiment::details { class Aging2Worker; }

namespace experiment::details {

class Aging2Master {
    friend class Aging2Worker;

    const Aging2Experiment& m_parameters;
    cuckoohash_map<uint64_t, uint64_t> m_vertices_final; // list of vertices in the final graph
    cuckoohash_map<uint64_t, bool> m_vertices_present; // current list of vertices present in the library
    cuckoohash_map<graph::Edge, bool> m_edges_present; // which edges are already present
    const bool m_is_directed; // is the graph directed?
    std::vector<Aging2Worker*> m_workers; // pool of workers
    std::atomic<int64_t> m_num_operations_performed = 0; // current number of operations performed so far
    std::atomic<int> m_last_progress_reported = 0; // the last progress of the experiment, reported by any of the worker threads. E.g. 1%, 2%, 3%, so on.
    uint64_t* m_vertices2remove = nullptr; // array with the vertices to remove, as they are not in the final graph, in the last stage of the experiment

    // Cumulative frequencies of the vertices
public: struct RankFrequency { uint64_t m_vertex_id; uint64_t m_frequency; };
private: RankFrequency* m_vertices_freq = nullptr;
    uint64_t m_vertices_freq_sz = 0;
    common::StaticIndex<uint64_t> m_vertices_freq_index;
    constexpr static uint64_t m_vertices_freq_index_leaf_sz = 64; // number of values per entry in the index

    // report how long it took to perform 1x, 2x, 3x, ... updates w.r.t. to the loaded graph.
    std::chrono::steady_clock::time_point m_time_start; // when the computation started
    uint64_t* m_reported_times = nullptr; // microsecs
    std::atomic<int> m_last_time_reported = 0;

    Aging2Result m_results; // final results of the experiment

    // Initialise the set of workers
    void init_workers();

    // Compute the cumulative frequencies for the vertices in the final graph. Here by frequency we mean the number of
    // outgoing edges attached to each node. The method may also insert nodes that do not exist in the final graph but
    // are basically additional `noise', according to the vertex expansion factor (m_parameters.m_ef_vertices)
    void init_cumulative_frequencies();

    // Execute the main part of the experiment, that is the insertions/deletions in the graph with the worker threads
    void do_run_experiment();

    // Remove the vertices that do not belong to the final graph
    void do_remove_artificial_vertices();

    // Save the current results in `m_results'
    void store_results();

    // print to stdout the number of vertices/edges expected and effectively stored in the final library. The two values should be equal.
    void log_num_vtx_edges();

public:
    Aging2Master(const Aging2Experiment& parameters);

    // Destructor
    ~Aging2Master();

    // Execute the actual experiment
    Aging2Result execute();

    // Is the graph directed?
    bool is_directed() const { return m_is_directed; }

    // Total number of vertices that can be generated
    uint64_t num_vertices() const { return m_vertices_freq_sz; };

    // Total number of operations to perform (insertions/deletions)
    uint64_t num_operations_total() const;

    // Max size, in terms of number of edges, that the intermediate graph can reach during the update phase
    uint64_t max_num_edges() const;

    // Return the maximum value in the cumulative sum of frequencies. That is the frequencies span in the range [0, distr_max_value())
    uint64_t distr_max_value() const;

    // Return the offset for the vertex id associated to the given cumulative frequency
    uint64_t distr_vertex_offset(uint64_t cumulative_frequency) const;

    // Access the configuration of this experiment
    const Aging2Experiment& parameters() const { return m_parameters; }
};

}

