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

#include "timer.hpp"

#include "time.hpp" // common::time::to_string()

using namespace std;
using namespace std::chrono;
using namespace common;

template <bool with_barrier>
string Timer<with_barrier>::to_string() const{
    // start and stop points in time
    if(m_t0 == clock::time_point{}) ERROR("Timer not even started");
    clock::time_point t1;
    if(m_t1 == clock::time_point{})
        t1 = clock::now();
    else
        t1 = m_t1;

    // duration
    auto d = t1 - m_t0;

    return time::to_string( d );
}

template<bool use_barrier>
std::ostream& (common::operator<<)(std::ostream& out, const common::Timer<use_barrier>& timer){
    out << timer.to_string();
    return out;
}

template<bool B1, bool B2>
Timer<B1> (common::operator+)(Timer<B1> timer1, Timer<B2> timer2){
    using clock = typename Timer<B1>::clock;
    typename clock::time_point nullpoint{};
    Timer<B1> result;

    // stop the timer if they are executing. Note we are passing the timers by value, no changes are forwarded to the original timers
    if(timer1.m_t0 != nullpoint && timer1.m_t1 == nullpoint){
        timer1.m_t1 = clock::now();
    }
    if(timer2.m_t0 != nullpoint && timer2.m_t1 == nullpoint){
        timer2.m_t1 = clock::now();
    }

    // first timer never launched?
    if(timer1.m_t0 == nullpoint){
        result.m_t0 = timer2.m_t0;
        result.m_t1 = timer2.m_t1;
    } else {
        // sum the duration from the two timers
        result = timer1;

        if(timer2.m_t0 != nullpoint){
            result.m_t0 -= (timer2.m_t1 - timer2.m_t0);
        }
    }

    return result;
};


// Template instantiantions
template class common::Timer<true>;
template class common::Timer<false>;
template ostream& (common::operator<<)(ostream& out, const Timer<true>& timer);
template ostream& (common::operator<<)(ostream& out, const Timer<false>& timer);
template Timer<true> (common::operator+)(Timer<true> t1, Timer<true> t2);
template Timer<true> (common::operator+)(Timer<true> t1, Timer<false> t2);
template Timer<false> (common::operator+)(Timer<false> t1, Timer<true> t2);
template Timer<false> (common::operator+)(Timer<false> t1, Timer<false> t2);
