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
#ifndef COMMON_SORTING_HPP
#define COMMON_SORTING_HPP

#include <cstdint>
#include <functional> // std::less

#include "details/sorting_impl.hpp"

namespace common {

    /**
     * Sort in parallel the input array T. The sorting algorithm is parallel.
     *
     * @param array the input array to sort
     * @param array_sz the size of the input array
     * @param fn_less a comparator function that returns true if the a < b.
     */
    template<typename T, typename FunctionLess>
    void sort(T* array, uint64_t array_sz, const FunctionLess& fn_less = std::less{}){
        details::sorting::implementation(array, array_sz, fn_less);
    }
} // namespace

#endif //COMMON_SORTING_HPP
