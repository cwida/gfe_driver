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

#include "timeout_service.hpp"

#include <cassert>

using namespace std;

namespace utility {

/*****************************************************************************
 *                                                                           *
 *  TimeoutService                                                           *
 *                                                                           *
 *****************************************************************************/
void TimeoutService::start() {
    assert(!m_background_thread.joinable() && "A background thread is already present");
    m_background_thread = thread(&TimeoutService::main_thread, this);
}

void TimeoutService::stop(){
    {
        scoped_lock<mutex> lock(m_mutex);
        m_terminate = true;
    }
    m_condvar.notify_all();

    m_background_thread.join(); // wait for termination
}

void TimeoutService::main_thread(){
    bool terminate = false;
    while(!terminate){
        unique_lock<mutex> lock(m_mutex);
        m_condvar.wait_for(lock, 1s, [this](){ return !m_terminate; });
        terminate = m_terminate;
        lock.unlock();

        m_is_timeout = ( clock_t::now() - m_start ) > m_budget;
    }
}

} // namespace
