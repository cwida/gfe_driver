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

#ifndef COMMON_OPTIMISATION_HPP
#define COMMON_OPTIMISATION_HPP

#include <cinttypes>

namespace common {

/**
 * Compiler barrier: at compiler level, avoid reordering instructions before and
 * after the barrier
 */
inline void compiler_barrier(){
    __asm__ __volatile__("": : :"memory");
};

/**
 * Read the CPU timestamp counter
 */
inline uint64_t rdtscp(){
    uint64_t rax;
    asm volatile (
    "rdtscp ; shl $32, %%rdx; or %%rdx, %%rax; "
    : "=a" (rax)
    : /* no inputs */
    : "rcx", "rdx"
    );
    return rax;
}

/**
 * Branch prediction macros
 */
#define LIKELY(x) __builtin_expect((x), 1)
#define UNLIKELY(x) __builtin_expect((x), 0)

}


#endif //COMMON_OPTIMISATION_HPP
