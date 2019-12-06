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

#include "build_thread.hpp"

#include "common/error.hpp"
#include "common/quantity.hpp" // for debugging purposes
#include "common/system.hpp"
#include "library/interface.hpp"

using namespace std;

/*****************************************************************************
 *                                                                           *
 * Debug                                                                     *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG
namespace gfe { extern mutex _log_mutex [[maybe_unused]]; }
#define COUT_DEBUG_FORCE(msg) { std::scoped_lock<std::mutex> lock{::gfe::_log_mutex}; std::cout << "[BuildThread::" << __FUNCTION__ << "] " << msg << std::endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif

/*****************************************************************************
 *                                                                           *
 * BuildThread impl                                                          *
 *                                                                           *
 *****************************************************************************/
namespace gfe::experiment::details {

BuildThread::BuildThread(std::shared_ptr<gfe::library::UpdateInterface> interface, int thread_id, std::chrono::milliseconds frequency) :
    m_interface(interface), m_thread_id(thread_id), m_frequency(frequency){
    m_terminate = true; // reset by the background thread
    if(m_frequency > 0ms){ // otherwise, never invoke #build()
        start();
    }
}

void BuildThread::start(){
    COUT_DEBUG("waiting for the service to start...");
    m_thread = thread{ &BuildThread::main_thread, this };

    unique_lock<mutex> lock(m_mutex);
    m_condvar.wait(lock, [this](){ return !m_terminate; });
    COUT_DEBUG("ack started");
}


BuildThread::~BuildThread(){
    stop();
}

void BuildThread::stop(){
    unique_lock<mutex> lock(m_mutex);
    if(!m_terminate){
        COUT_DEBUG("waiting for the service to stop...");
        m_terminate = true;
        lock.unlock();
        m_condvar.notify_all();
        m_thread.join();
        m_condvar.notify_all();
        COUT_DEBUG("ack terminated");
    } else { // in case multiple threads invoked #stop
        m_condvar.wait(lock, [this](){ return !m_thread.joinable(); });
    }
}

void BuildThread::main_thread(){
    COUT_DEBUG("service started, thread_id: " << m_thread_id << ", frequency: " << common::DurationQuantity(m_frequency));
    common::concurrency::set_thread_name("build service");

    unique_lock<mutex> lock(m_mutex);
    bool terminate = m_terminate = false;
    lock.unlock();
    m_condvar.notify_all();

    m_interface->on_thread_init(m_thread_id);

    do {
        lock.lock();
        m_condvar.wait_for(lock, m_frequency, [this](){ return m_terminate; });
        terminate = m_terminate;
        lock.unlock();

        // no need to hold the lock here
        COUT_DEBUG("#build, num invocations: " << m_num_invocations << ", terminate: " << boolalpha << terminate);
        m_interface->build();
        m_num_invocations++;
    } while(!terminate);

    m_interface->on_thread_destroy(m_thread_id);

    COUT_DEBUG("service terminated");
}

} // namespace
