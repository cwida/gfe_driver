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

#include "graph/edge.hpp"
#include "library/interface.hpp"

#include <condition_variable>
#include <mutex>
#include <thread>

namespace gfe::experiment::details {

/**
 * This class is not thread safe.
 */
class AsyncBatch {
    gfe::library::UpdateInterface* m_interface;
    gfe::library::UpdateInterface::SingleUpdate** m_batches = nullptr; // array of batches
    size_t* m_batches_num_entries = nullptr; // number of filled entries in each batch
    const int m_batches_sz; // number of batches
    const int m_batch_sz; // the size of each batch, as multiples of library::UpdateInterface::SingleUpdate
    int m_batch_index = 0; // the current batch being filled
    int m_batch_pos = 0; // the next position free in the current batch
    const int m_thread_id; // the thread id in the interface
    std::thread m_thread_handle; // the handle to the thread invoking the interface
    std::mutex m_thread_mutex; // sync between the producer and the consumer
    std::condition_variable m_thread_condvar;
    int m_send_index = 0; // the current batch being processed asynchronously
    int m_send_upto = 0; // the number of batches to process asynchronously, only read & update while holding `m_thread_mutex'

    // Asychronous thread
    void main_thread();

public:
    /**
     * Constructor
     * @param interface used to import the batches in the graph
     * @param thread_id the thread_id passed to the interface (on_worker_init)
     * @param num_batches the number of batches that can be asynchronously processed
     * @param batch_sz the size of each batch
     */
    AsyncBatch(gfe::library::UpdateInterface* interface, int thread_id, int num_batches, int batch_sz);

    /**
     * Destructor. It implicitly flushes any pending updates and stops the service.
     */
    ~AsyncBatch();

    /**
     * Add the given edge to the graph
     */
    void add_edge(gfe::graph::WeightedEdge edge);


    /**
     * Remove the given edge from the graph
     */
    void remove_edge(gfe::graph::Edge edge);

    /**
     * Prepare the next batch
     * @param synchronise if true, wait for all pending batches to be processed, before resuming
     */
    void flush(bool synchronise);

};

} // namespace


