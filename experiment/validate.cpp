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

#include "validate.hpp"

#include "graph/edge_stream.hpp"
#include "library/interface.hpp"
#include "configuration.hpp"

using namespace std;

namespace gfe::experiment {

uint64_t validate_updates(shared_ptr<gfe::library::Interface> ptr_interface, shared_ptr<gfe::graph::WeightedEdgeStream> ptr_stream) {
    auto interface = ptr_interface.get();
    auto stream = ptr_stream;

    LOG("Validation started");

    uint64_t num_threads = thread::hardware_concurrency();
    interface->on_main_init(num_threads);
    atomic<int64_t> num_errors = 0;

    auto routine = [stream, interface, &num_errors](int thread_id, uint64_t from, uint64_t to){
        interface->on_thread_init(thread_id);

        for(uint64_t i = from; i < to; i++){
            auto edge = stream->get(i);
            auto w1 = interface->get_weight(edge.source(), edge.destination());
            if(w1 != edge.m_weight){
                LOG("ERROR [" << i << "] Edge mismatch " << edge.source() << " -> " << edge.destination() << ", retrieved weight: " << w1 << ", expected: " << edge.weight());
                num_errors++;
            }
            if(interface->is_undirected()){
                auto w2 = interface->get_weight(edge.destination(), edge.source());
                if(w2 != edge.m_weight){
                    LOG("ERROR [" << i << "] Edge mismatch " << edge.source() << " <- " << edge.destination() << ", retrieved weight: " << w1 << ", expected: " << edge.weight());
                    num_errors++;
                }
            }
        }
        interface->on_thread_destroy(thread_id);
    };


    uint64_t edges_per_thread = stream->num_edges() / num_threads;
    uint64_t odd_threads = stream->num_edges() % num_threads;

    vector<thread> threads;
    uint64_t from = 0;
    for(uint64_t i = 0; i < num_threads; i++){
        uint64_t to = from + edges_per_thread + (i < odd_threads);
        threads.emplace_back(routine, i, from, to);
        from = to;
    }

    for(auto& t: threads) t.join();

    interface->on_main_destroy();

    if(num_errors == 0){
        LOG("Validation succeeded");
    } else {
        LOG("Number of validation errors: " << num_errors);
    }

    return num_errors;
}

} // namespace



