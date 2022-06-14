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

#ifndef COMMON_PERMUTATION_HPP
#define COMMON_PERMUTATION_HPP

#include <cinttypes>

namespace common {

/**
 * Create a random permutation of the allocated (but not initialised) array of cardinality array_sz.
 * The final elements will be a shuffle of the `array_sz' elements in [0, array_sz).
 */
template<typename T>
void permute(T* array, uint64_t array_sz, uint64_t seed);

}

#endif //COMMON_PERMUTATION_HPP
