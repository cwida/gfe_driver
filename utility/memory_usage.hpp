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

#include <cstdint>
#include <string>
#include <vector>

namespace gfe::utility {

/**
 * Compute the memory usage of this process.
 */
class MemoryUsage {
public:
    /**
     * This method must be invoked at start-up. The process is reloaded with the functions
     * #malloc/#free override by the loader.
     * The function can be invoked multiple times, subsequent calls to the method are transformed into *nop*.
     *
     * This method is not thread-safe.
     */
    static void initialise(int argc, char* argv[]);

    /**
     * Check whether the instrumentation facility has been initialised.
     */
    static bool is_initialised();

    /**
     * Compute the memory footprint of this process.
     * For allocations made with #malloc & co, by the C++ runtime, it reports the virtual memory allocated.
     * For allocations explicitly made with #mmap, it reports the physical memory currently in use.
     *
     * The function always returns `0` if the class has not been previously initialised.
     */
    static int64_t memory_footprint();

    /**
     * Get the virtual space used by the given allocation, assuming the allocation has been made by glibc
     */
    static uint64_t get_allocated_space(const void* pointer);
};

} // namespace
