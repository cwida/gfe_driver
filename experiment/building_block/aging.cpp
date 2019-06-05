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

#include "aging.hpp"

#include <cassert>
#include <random>
#include <thread>
#include <vector>

#include "common/optimisation.hpp"
#include "common/timer.hpp"
#include "graph/edge.hpp"
#include "graph/edge_stream.hpp"
#include "library/interface.hpp"
#include "configuration.hpp"

using namespace common;
using namespace std;

namespace experiment {

Aging::Aging(std::shared_ptr<library::UpdateInterface> interface) : Aging(interface, make_shared<graph::WeightedEdgeStream>()){ }
Aging::Aging(std::shared_ptr<library::UpdateInterface> interface, std::shared_ptr<graph::WeightedEdgeStream> stream) : Aging(interface, stream, configuration().num_threads(THREADS_WRITE)) { }
Aging::Aging(std::shared_ptr<library::UpdateInterface> interface, std::shared_ptr<graph::WeightedEdgeStream> stream, int64_t num_threads) : m_interface(interface), m_stream(stream), m_num_threads(num_threads) { }

void Aging::set_expansion_factor(double factor){
    if(factor < 1) INVALID_ARGUMENT("the expansion factor must be >= 1, instead the value given is: " << factor);
    m_expansion_factor = factor;
}

void Aging::set_operation_granularity(uint64_t granularity){
    if(granularity < 1) INVALID_ARGUMENT("the granularity given must be > 0: " << granularity);
    m_granularity = granularity;
}

void Aging::thread_main(int thread_id){
    assert(thread_id >= 0 && "The thread id should be >= 0");
    auto interface = m_interface.get();
    interface->on_thread_init(thread_id);

    // internal generator
    mt19937_64 random( random_device() );
    uniform_real_distribution<double> uniform{ 0., 1. };
    vector<graph::Edge> edges2remove;

    // define the range of this partition
    const uint64_t gen_part_len = m_stream->num_edges() / m_num_threads;
    const uint64_t gen_part_mod = m_stream->num_edges() % m_num_threads;
    const uint64_t part_start = thread_id * gen_part_len + min<int>(gen_part_mod, thread_id); // incl.
    const uint64_t part_length = gen_part_len + (thread_id < gen_part_mod);
    const uint64_t part_end = part_start + part_length; // excl.
    uint64_t part_pos = part_start; // current position in the range position

    // wait for all threads to init
    compiler_barrier();
    m_startup_counter--;
    while(m_startup_counter > 0) /* nop */;
    compiler_barrier();


    // current count of operations performed
    int64_t global_operation_count = 0;
    while( (global_operation_count = m_num_operations_performed.fetch_add(m_granularity)) < m_num_operations_total ){
        // shall we perform a burst of insertions or deletions ?
        if(interface->num_edges() < static_cast<uint64_t>(m_expansion_factor * m_stream->num_edges()) &&
                (m_num_operations_total - global_operation_count) < static_cast<int64_t>(edges2remove.size())){
            // perform `m_granularity' insertions ...

            for(uint64_t i = 0; i < m_granularity; i++){
                assert(part_pos <= part_end);
                uint64_t missing_edges_from_final_graph = part_end - part_pos;
                if ( uniform(random) < static_cast<double>(missing_edges_from_final_graph) / (missing_edges_from_final_graph + (m_num_operations_total - global_operation_count)) ) {
                    assert(part_pos < part_end);

                    auto weighted_edge = m_stream->get(part_pos);
                    auto edge = weighted_edge.edge();

                    if( m_edges_present.insert_or_assign(edge, FINAL ) ) {
                        thread_insert(weighted_edge);
                    }



                    part_pos++;
                }


                global_operation_count++;
            }



        }

    }




    interface->on_thread_destroy(thread_id);
}



// run the experiment
std::chrono::microseconds Aging::execute(){
    auto interface = m_interface.get();
    // start the threads
    LOG("Initialising " << m_num_threads << " threads ...");
    vector<thread> threads;
    for(int64_t i = 0; i < m_num_threads; i++){
        threads.emplace_back(&Aging::thread_main, this, static_cast<int>(i));
    }

    // wait for all threads to init
    LOG("Waiting for all threads to start...");
    compiler_barrier();
    while(m_startup_counter > 0) /* nop */;
    compiler_barrier();

    // start the timer
    Timer timer;
    timer.start();
    LOG("Experiment started!");


    // wait for all threads to finish
    for(auto& t : threads) t.join();
    timer.stop();

    LOG("Experiment terminated in " << timer);


    return timer.duration<chrono::microseconds>();
}




} // namespace experiment
