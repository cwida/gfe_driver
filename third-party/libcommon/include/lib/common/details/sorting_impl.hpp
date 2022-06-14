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
#ifndef COMMON_SORTING_IMPL_HPP
#define COMMON_SORTING_IMPL_HPP

#include <algorithm>
#include <cstdint>
#include <future>
#include <memory>
#include <thread>
#include <vector>

#include "../sampling.hpp"

namespace common::details::sorting {

    // Partition the `array' according to the given samples
    template<typename T, typename FunctionLess>
    void partition(T*  __restrict array, uint64_t array_sz, const FunctionLess& fn_less, const T* __restrict samples, uint64_t* __restrict intervals, uint64_t num_samples) {
        for(uint64_t i = 0; i < array_sz; i++){
            int64_t k = num_samples -1;
            uint64_t position = i;
            while(k >= 0 && fn_less(array[position], samples[k]) ){
                if(position > intervals[k]){
                    std::swap(array[position], array[intervals[k]]);
                    position = intervals[k];
                }

                intervals[k]++;

                // next iteration
                k--;
            }
        }
    }

    // Sort the input array
    // Based on the sample sort procedure of Section § 5.7.2 in
    // K. Mehlhorn, P. Sanders, Algorithms and Data Structures. The Basic Toolbox, Springer 2008.
    template<typename T, typename FunctionLess>
    void implementation(T* array, uint64_t array_sz, const FunctionLess& fn_less, uint64_t num_samples = std::thread::hardware_concurrency()){
        // samples => O(k), k = num samples
        if(num_samples >= array_sz){ std::sort(array, array + array_sz, fn_less); return; }
        std::unique_ptr<T[]> ptr_samples { new T[num_samples]() };
        std::unique_ptr<uint64_t[]> ptr_intervals { new uint64_t[num_samples]() };
        T* samples = ptr_samples.get();
        uint64_t* intervals = ptr_intervals.get();
        common::random_sample(array, array_sz, samples, num_samples);
        std::sort(samples, samples + num_samples, fn_less);

        // sequentially partition the intervals, this is the heavy part of the algorithm => O(n), n = array_sz
        partition(array, array_sz, fn_less, samples, intervals, num_samples);

        // sort the intervals independently in parallel
        std::vector<std::future<void>> tasks;
        uint64_t start = 0; // inclusive
        for(uint64_t i = 0; i < num_samples; i++){
            uint64_t end = intervals[i]; // exclusive
            if(start < end){ // if the interval is not empty
                tasks.push_back(std::async(std::launch::async, [&array, start, end, &fn_less](){
                    std::sort(array + start, array + end, fn_less);
                }));
            }
            start = end; // next iteration
        }
        if(start < array_sz){ // last interval
            uint64_t end = array_sz;
            tasks.push_back(std::async(std::launch::async, [&array, start, end, &fn_less](){
                std::sort(array + start, array + end, fn_less);
            }));
        }

        for(auto& t: tasks){ t.get(); }; // wait for all the tasks to complete
    }

} // namespace

#endif //COMMON_SORTING_IMPL_HPP
