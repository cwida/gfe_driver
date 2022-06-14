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

#ifndef COMMON_MATH_HPP
#define COMMON_MATH_HPP

#include <cinttypes>

namespace common {

/**
 * Get the lowest power of 2, p,  such that p >= x
 * That is: 2^ceil(log2(x))
 */
uint64_t hyperceil(uint64_t x);


} // common

#endif //COMMON_MATH_HPP
