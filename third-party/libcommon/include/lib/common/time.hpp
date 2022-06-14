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

#ifndef COMMON_TIME_HPP
#define COMMON_TIME_HPP

#include <cinttypes>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>

namespace common::time {

template <typename D>
std::string to_nanoseconds(D duration){
    using namespace std::chrono;
    std::stringstream result;
    result << (uint64_t) duration_cast<std::chrono::nanoseconds>(duration).count() << " ns";
    return result.str();
}

template <typename D>
std::string to_microseconds(D duration){
    using namespace std::chrono;
    uint64_t time_in_nanosecs = (uint64_t) duration_cast<std::chrono::nanoseconds>(duration).count();
    uint64_t time_in_microsecs = time_in_nanosecs / 1000;

    std::stringstream result;
    if(time_in_microsecs >= 3){
        result << time_in_microsecs << " us";
    } else {
        char buffer[128];
        snprintf(buffer, 128, "%.3d", (int) (time_in_nanosecs % 1000));
        result << time_in_microsecs << "." << buffer << " us";
    }

    return result.str();
}

template <typename D>
std::string to_milliseconds(D duration){
    using namespace std::chrono;
    uint64_t time_in_microsecs = (uint64_t) duration_cast<std::chrono::microseconds>(duration).count();
    uint64_t time_in_millisecs = time_in_microsecs / 1000;

    std::stringstream result;
    if(time_in_microsecs >= 3){
        result << time_in_millisecs << " ms";
    } else {
        char buffer[128];
        snprintf(buffer, 128, "%.3d", (int) (time_in_microsecs % 1000));
        result << time_in_millisecs << "." << buffer << " ms";
    }

    return result.str();
}

template <typename D>
std::string to_seconds(D duration){
    using namespace std::chrono;
    uint64_t time_in_millisecs = (uint64_t) duration_cast<std::chrono::milliseconds>(duration).count();
    uint64_t time_in_seconds = time_in_millisecs / 1000;

    std::stringstream result;
    char buffer[128];
    snprintf(buffer, 128, "%.3d", (int) (time_in_millisecs % 1000));
    result << time_in_seconds << "." << buffer << " s";

    return result.str();
}

template <typename D>
std::string to_minutes(D duration){
    using namespace std::chrono;
    uint64_t seconds = ((uint64_t) duration_cast<std::chrono::seconds>(duration).count()) % 60ull;
    uint64_t minutes = (uint64_t) duration_cast<std::chrono::minutes>(duration).count();

    char buffer[128];
    snprintf(buffer, 128, "%" PRIu64 ":%02" PRIu64 " mins", minutes, seconds);
    return std::string(buffer);
}

template <typename D>
std::string to_hours(D duration){
    using namespace std::chrono;
    uint64_t seconds = ((uint64_t) duration_cast<std::chrono::seconds>(duration).count()) % 60ull;
    uint64_t minutes = (uint64_t) duration_cast<std::chrono::minutes>(duration).count() % 60ull;
    uint64_t hours = (uint64_t) duration_cast<std::chrono::hours>(duration).count();

    char buffer[128];
    snprintf(buffer, 128, "%" PRIu64 ":%02" PRIu64 ":%02" PRIu64 " hours", hours, minutes, seconds);
    return std::string(buffer);
}

template <typename D>
std::string to_days(D duration){
    using namespace std::chrono;
    uint64_t seconds = ((uint64_t) duration_cast<std::chrono::seconds>(duration).count()) % 60ull;
    uint64_t minutes = (uint64_t) duration_cast<std::chrono::minutes>(duration).count() % 60ull;
    uint64_t hours = (uint64_t) duration_cast<std::chrono::hours>(duration).count() % 24ull;
    uint64_t days = (uint64_t) duration_cast<std::chrono::hours>(duration).count() / 24;

    char buffer[128];
    snprintf(buffer, 128, "%" PRIu64 " day(s) and %02" PRIu64 ":%02" PRIu64 ":%02" PRIu64 " hours", days, hours, minutes, seconds);
    return std::string(buffer);
}

template <typename Duration>
std::string to_string(Duration d){
    uint64_t time_in_nanosecs = (uint64_t) std::chrono::duration_cast<std::chrono::nanoseconds>(d).count();
    if(time_in_nanosecs <= 1000){
        return to_nanoseconds(d);
    } else if(time_in_nanosecs <= (uint64_t) std::pow(10, 6)){
        return to_microseconds(d);
    } else if(time_in_nanosecs <= (uint64_t) std::pow(10, 9)) {
        return to_milliseconds(d);
    } else if(time_in_nanosecs <= (uint64_t) std::pow(10, 9) * 90){ // 90 seconds
        return to_seconds(d);
    } else if(time_in_nanosecs < (uint64_t) std::pow(10, 9) * 60 * 60){
        return to_minutes(d);
    } else if(time_in_nanosecs < (uint64_t) std::pow(10, 9) * 60 * 60 * 24){
        return to_hours(d);
    } else {
        return to_days(d);
    }
}

} // namespace

#endif // COMMON_TIME_HPP
