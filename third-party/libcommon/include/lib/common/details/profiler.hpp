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

#ifndef COMMON_PROFILER_DETAILS_HPP
#define COMMON_PROFILER_DETAILS_HPP

#include <cstddef>

namespace common { namespace details {

/**
 * Initialise the PAPI library (if not already initialised) and retrieve the event code
 * associated to a named event
 */
class BaseProfiler {
    static bool library_initialised;
    static void initialise_library();

public:
    BaseProfiler();

    /**
     * Return -1 if the event is not available, otherwise it's PAPI event code
     */
    static int get_event_code(const char* event_name);
};


/**
 * Boiler plate to register the PAPI events and start and stop the recording of perf_events.
 */
class GenericProfiler : public BaseProfiler {
    static constexpr int m_events_capacity = 8;

protected:
    int m_events[m_events_capacity];
    int m_events_sz = 0;
    int m_event_set = 0;

    void add_events(const char* errorstring, const char* event_name);
    void add_events(const char* errorstring, const char* alternative_events[], size_t num_alternative_events);
    void register_events();
    void unregister_events();

    void start();
    void stop(long long* resultset);
    void snapshot(long long* resultset);

public:
    GenericProfiler();
    ~GenericProfiler();
};

}} // common::details

#endif //COMMON_PROFILER_DETAILS_HPP
