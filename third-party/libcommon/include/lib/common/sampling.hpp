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

#ifndef COMMON_SAMPLING_HPP
#define COMMON_SAMPLING_HPP

#include <cinttypes>
#include <random>

#include "details/sampling_impl.hpp"

namespace common {

/**
 * Retrieve a random sample of the input following a uniform distribution
 *
 * @param input the array to sample
 * @param input_sz the size of the array `input'
 * @param output the output array where the samples will be copied. It must be preallocated by the caller with a
 *        capacity of at least of `num_samples'
 * @param num_samples the total number of samples to retrieve
 * @param seed the seed for the random generator
 */
template<typename T>
void random_sample(const T* input, uint64_t input_sz, T* output, uint64_t& num_samples, uint64_t seed = std::random_device{}()){
    details::sampling::implementation(input, input_sz, output, num_samples, seed); // details/sampling_impl.hpp
}

} // namespace

#endif //COMMON_SAMPLING_HPP
