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

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <memory>
#include <mutex>
#include <thread>

namespace gfe::library { class UpdateInterface; } // forward decl.

namespace gfe::experiment::details {

/**
 * This service continuously invoke the library's method #build() every tot seconds.
 * The idea is, on delta-based systems, that every tot seconds a new snapshot is built
 * every tot seconds.
 * On all the other systems, an invocation to #build() becomes a nop.
 */
class BuildThread {
    BuildThread(const BuildThread&) = delete;
    BuildThread& operator=(const BuildThread&) = delete;

    std::shared_ptr<gfe::library::UpdateInterface> m_interface; // the library where to invoke the method #build
    const int m_thread_id; // the internal thread_id to use with #on_thread_init and #on_thread_exit
    const std::chrono::milliseconds m_frequency; // how frequently the service shall invoke #build()
    std::atomic<uint64_t> m_num_invocations = 0; // the total number of calls to #build() by the service, so far

    bool m_terminate = false; // signal the background thread that 1) has started and 2) has terminated
    std::mutex m_mutex; // sync to start/terminate the service
    std::condition_variable m_condvar; // as above
    std::thread m_thread; // the handle for the thread used by the background service

    // start the service, that is the background thread
    void start();

    // the actual logic of the background thread
    void main_thread();

public:
    /**
     * Constructor. It implicitly starts the service/background thread invoking #build
     * @param interface the library where to invoke the method #build
     * @param thread_id the thread_id passed to the library and used by the service/background thread
     * @param frequency how frequently the method #build() shall be invoked
     */
    BuildThread(std::shared_ptr<gfe::library::UpdateInterface> interface, int thread_id, std::chrono::milliseconds frequency);

    /**
     * Destructor. It implicitly stops the service.
     */
    ~BuildThread();

    /**
     * It stops the service, that is, the background thread.
     */
    void stop();

    /**
     * Retrieve the total number of invocations to #build() so far
     */
    uint64_t num_invocations() const { return m_num_invocations; }
};

} // namespace

