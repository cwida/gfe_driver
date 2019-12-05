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

#include <chrono>
#include <condition_variable>
#include <mutex>
#include <thread>

namespace gfe::utility {
    
/**
 * This service keeps track sets the flag is_timeout() after a certain given of time has passed
 * since the service itself was created.
 * The service is meant to be used by multiple threads to poll continuously whether they can
 * continue their computation or they depleted their budget and abort the computation.
 * This class is thread safe.
 */
class TimeoutService {
    using clock_t = std::chrono::steady_clock;
    const clock_t::time_point m_start; // the time when the service was started
    const std::chrono::seconds m_budget; // the amount of time that must pass before updating the flag `m_is_timeout'
    bool m_is_timeout = false; // the flag to update asynchronously
    bool m_terminate = false; // signal termination to the async thread
    std::mutex m_mutex; // sync the termination with the background thread
    std::condition_variable m_condvar; // sync the termination with the background thread
    std::thread m_background_thread; // handle to the background thread

    // The underlying thread responsible to update asynchronously the flag `m_is_timeout'
    void main_thread();

    // Starts the service
    void start();

    // Stops the service
    void stop();
public:
    /**
     * Creates the service.
     * Note: giving the value timeout = 0 has the special behaviour that the service will never start,
     * and the method #is_timeout() will always return false. It's like an indefinite time budget.
     *
     * @param timeout the amount of time that must pass before the method is_timeout() can return true;
     */
    TimeoutService(std::chrono::seconds timeout) : m_start(clock_t::now()), m_budget(timeout){
        start();
    };

    /**
     * Destructor. It implicitly stops the service
     */
    ~TimeoutService(){ stop(); };

    /**
     * Checks whether the specified amount of time has passed.
     */
    bool is_timeout() const { return m_is_timeout; }
};
    
} // namespace
