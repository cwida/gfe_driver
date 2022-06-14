/**
 * This file is part of libcommon.
 *
 * libcommon is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, orF
 * (at your option) any later version.
 *
 * libcommon is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with libcommon.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef COMMON_TIMER_HPP
#define COMMON_TIMER_HPP

#include <chrono>
#include <ostream>
#include "error.hpp"
#include "optimisation.hpp"

namespace common {

template<bool use_barrier = true>
class Timer {
    template<bool B1, bool B2> friend Timer<B1> operator+(Timer<B1>, Timer<B2>);
//    Timer(const Timer&) = delete;
//    Timer& operator=(const Timer& timer) = delete;
    using clock = std::chrono::steady_clock;

    clock::time_point m_t0; // start time
    clock::time_point m_t1; // end time

public:
    Timer(){ }

    void start(){
        m_t1 = clock::time_point{};
        if(use_barrier) compiler_barrier();
        m_t0 = clock::now();
        if(use_barrier) compiler_barrier();
    }

    void resume(){
        if(m_t0 != clock::time_point{}) return; // already running;
        if(m_t1 == clock::time_point{}){ // this timer has never been executed
            start();
        } else {
            if(use_barrier) compiler_barrier();
            m_t0 = clock::now() - (m_t1 - m_t0);
            if(use_barrier) compiler_barrier();
        }
    }

    void stop() {
        if(use_barrier) compiler_barrier();
        m_t1 = clock::now();
        if(use_barrier) compiler_barrier();
    }

    template<typename D>
    D duration() const {
        return std::chrono::duration_cast<D>(m_t1 - m_t0);
    }

    template<typename D>
    uint64_t convert() const {
        return static_cast<uint64_t>( duration<D>().count() );
    }

    uint64_t nanoseconds() const{ return convert<std::chrono::nanoseconds>(); }

    uint64_t microseconds() const{ return convert<std::chrono::microseconds>(); }

    uint64_t milliseconds() const{ return convert<std::chrono::milliseconds>(); }

    uint64_t seconds() const{ return convert<std::chrono::seconds>(); }

    std::string to_string() const;

};

template<bool use_barrier>
std::ostream& operator<<(std::ostream& out, const Timer<use_barrier>& timer);

template<bool t1_barrier, bool t2_barrier>
common::Timer<t1_barrier> operator+(common::Timer<t1_barrier> t1, common::Timer<t2_barrier> t2);

} // namespace common

#endif //COMMON_TIMER_HPP
