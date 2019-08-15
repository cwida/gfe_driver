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

#include "async_batch.hpp"

#include <cassert>

#include "configuration.hpp" // LOG
#include "common/system.hpp"

using namespace std;

namespace experiment::details {

/*****************************************************************************
 *                                                                           *
 * Debug                                                                     *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG
#define COUT_DEBUG_FORCE(msg) { LOG("[AsyncBatch::" << __FUNCTION__ << "] [" << common::concurrency::get_thread_id() << "] " << msg); }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif

/*****************************************************************************
 *                                                                           *
 * Constructor                                                               *
 *                                                                           *
 *****************************************************************************/

AsyncBatch::AsyncBatch(library::UpdateInterface* interface, int thread_id, int num_batches, int batch_sz) : m_interface(interface), m_batches_sz(num_batches), m_batch_sz(batch_sz), m_thread_id(thread_id) {
    if(num_batches < 1){ throw std::invalid_argument("num_batches < 1"); }
    if(batch_sz < 1){ throw std::invalid_argument("batch_sz < 1"); }

    m_batches = new library::UpdateInterface::SingleUpdate*[m_batches_sz];
    m_batches_num_entries = new size_t[m_batch_sz]();
//    m_send_futures = new std::future<void>*[m_batch_sz](); // init the entries to nullptr
    for(int batch_index = 0; batch_index < m_batches_sz; batch_index++){
        m_batches[batch_index] = new library::UpdateInterface::SingleUpdate[m_batch_sz];
    }

    // start the thread
    unique_lock<mutex> lock(m_thread_mutex);
    m_send_upto = 1; // use to sync with the thread
    m_thread_handle = std::thread(&AsyncBatch::main_thread, this);
    m_thread_condvar.wait(lock, [this](){ return m_send_upto == 0; });
}

AsyncBatch::~AsyncBatch() {
    flush(true);

    { // wait for the async thread to terminate
        unique_lock<mutex> lock(m_thread_mutex);
        assert(m_send_upto == 0); // because of flush(true);
        m_send_upto = -1;
    }
    m_thread_condvar.notify_all();
    m_thread_handle.join();

    for(int batch_index = 0; batch_index < m_batches_sz; batch_index++){
        delete[] m_batches[batch_index]; m_batches[batch_index] = nullptr;
    }

    delete[] m_batches; m_batches = nullptr;
    delete[] m_batches_num_entries; m_batches_num_entries = nullptr;
//    delete[] m_send_futures; m_send_futures = nullptr;
}


/*****************************************************************************
 *                                                                           *
 * Async thread                                                              *
 *                                                                           *
 *****************************************************************************/

void AsyncBatch::main_thread(){
    m_interface->on_thread_init(m_thread_id);

    { // wake up the ctor that started this thread
        unique_lock<mutex> lock(m_thread_mutex);
        m_send_upto = 0;
    }
    m_thread_condvar.notify_all();
    COUT_DEBUG("Thread initialised");

    bool terminate = false;
    { // wait for a job to be available
        unique_lock<mutex> lock(m_thread_mutex);
        m_thread_condvar.wait(lock, [this](){ return m_send_upto != 0; });
        if(m_send_upto == -1){ // terminate
            m_send_upto = -2;
            terminate = true;
        } else {
            assert(m_send_upto > 0);
        }
    }

    while(!terminate){
        COUT_DEBUG("Process batch at index " << m_send_index << ", size: " << m_batches_num_entries[m_send_index]);
        m_interface->batch(m_batches[m_send_index], m_batches_num_entries[m_send_index], /* force */ true);
        m_send_index = (m_send_index +1) % m_batches_sz; // next batch to send
        m_thread_condvar.notify_all();

        { // wait for the next job
            unique_lock<mutex> lock(m_thread_mutex);
            m_send_upto--; // we just processed a job
            m_thread_condvar.wait(lock, [this](){ return m_send_upto != 0; });
            if(m_send_upto == -1){ // terminate
                m_send_upto = -2;
                terminate = true;
            } else {
                assert(m_send_upto > 0);
            }
        }
    }

    m_thread_condvar.notify_all(); // we are done
    m_interface->on_thread_destroy(m_thread_id);
    COUT_DEBUG("Thread terminated");
}


void AsyncBatch::flush(bool synchronise){
    unique_lock<mutex> lock(m_thread_mutex);
    COUT_DEBUG("sync: " << synchronise << ", m_send_up: " << m_send_upto);

    if(m_send_upto >= (m_batches_sz -1)){ // already full, wait for one batch to be done
        m_thread_condvar.wait(lock, [this](){ return m_send_upto < m_batches_sz -1; });
    }

    if(m_batch_pos > 0){ // send the next batch
        m_batches_num_entries[m_batch_index] = m_batch_pos;
        m_batch_pos = 0;
        m_batch_index = (m_batch_index + 1) % m_batches_sz;
        m_send_upto++;
    }

    lock.unlock();
    m_thread_condvar.notify_all(); // tell the async thread there are updates to process

    if(synchronise){ // wait for all batches to be sent
        lock.lock();
        m_thread_condvar.wait(lock, [this](){ return m_send_upto == 0; });
    }
}


/*****************************************************************************
 *                                                                           *
 * Edge insert/remove                                                        *
 *                                                                           *
 *****************************************************************************/

void AsyncBatch::add_edge(graph::WeightedEdge edge) {
    if(m_batch_pos == m_batch_sz) flush(false); // send all pending updates in the batch, reset m_batch_pos to 0
    assert(m_batch_pos < m_batch_sz && "Overflow");
    COUT_DEBUG("add_edge: " << edge << ", batch_index: " << m_batch_index << ", batch_pos: " << m_batch_pos);

    library::UpdateInterface::SingleUpdate* __restrict update = m_batches[m_batch_index] + m_batch_pos;
    update->m_source = edge.m_source;
    update->m_destination = edge.m_destination;
    update->m_weight = edge.m_weight;
    m_batch_pos++;
}

void AsyncBatch::remove_edge(graph::Edge edge) {
    if(m_batch_pos == m_batch_sz) flush(false); // send all pending updates in the batch, reset m_batch_pos to 0
    assert(m_batch_pos < m_batch_sz && "Overflow");

    library::UpdateInterface::SingleUpdate* __restrict update = m_batches[m_batch_index] + m_batch_pos;
    update->m_source = edge.m_source;
    update->m_destination = edge.m_destination;
    update->m_weight = -1;
    m_batch_pos++;
}


} // namespace


