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
 
#include "permutation.hpp"

#include <cassert>
#include <future>
#include <iostream>
#include <mutex>
#include <random>
#include <thread>
#include <vector>

using namespace std;

// This implementation is from the C++ driver for the experiments with the packed memory arrays.
// Which I believe is in turn derived from a paper from P. Sanders, needing to be proper referenced.
namespace common {

namespace {

template<typename T>
struct Bucket{
    using random_generator_t = mt19937_64;

    const int m_bucket_no;
    random_generator_t m_random_generator;
    vector<T>** m_chunks;
    const int m_chunks_size;
    uint64_t m_permutation_sz; // the size of the local permutation
    T* m_permutation;

    Bucket(int id, uint64_t seed, uint64_t no_chucks) :
            m_bucket_no(id), m_random_generator(seed + m_bucket_no), m_chunks(nullptr), m_chunks_size(no_chucks),
            m_permutation_sz(0), m_permutation(nullptr){
        // initialise the chunks
        m_chunks = new vector<T>*[m_chunks_size];
        for(int i = 0; i < m_chunks_size; i++){
            m_chunks[i] = new vector<T>();
        }
    }

    ~Bucket(){
        if(m_chunks != nullptr){
            for(int i = 0; i < m_chunks_size; i++){
                delete m_chunks[i]; m_chunks[i] = nullptr;
            }
        }
        delete[] m_chunks; m_chunks = nullptr;

        // DO NOT DELETE m_permutation, THE CLASS DOES NOT OWN THIS PTR
        m_permutation = nullptr;
    }
};

} // anonymous namespace

template<typename T>
static void do_permute(T* array, uint64_t array_sz, uint64_t no_buckets, uint64_t seed){
    if(no_buckets > array_sz) no_buckets = array_sz;

    // initialise the buckets
    vector< unique_ptr< Bucket<T> > > buckets;
//    size_t bytes_per_element = compute_bytes_per_elements(size); // the original implementation used an array with variable-length data
    for(uint64_t i = 0; i < no_buckets; i++){
        buckets.emplace_back(new Bucket<T>(i, seed + i, no_buckets /*, bytes_per_element*/));
    }

    { // exchange the elements in the buckets
        auto create_partition = [&buckets](int bucket_id, uint64_t range_start /* inclusive */ , uint64_t range_end /* exclusive */){
            assert(bucket_id < (int) buckets.size());

            // first gather all the buckets
            vector< vector<T>* > stores;
            for(size_t i = 0; i < buckets.size(); i++){
                Bucket<T>* b = buckets[i].get();
                stores.push_back( b->m_chunks[bucket_id] );
            }

            // current state
            Bucket<T>* bucket = buckets[bucket_id].get();
            uniform_int_distribution<size_t> distribution{0, buckets.size() -1};
            for(uint64_t i = range_start; i < range_end; i++){
                size_t target_bucket = distribution(bucket->m_random_generator);
                stores[target_bucket]->push_back(i);
            }
        };

        std::vector<future<void>> tasks;
        size_t range_start = 0, range_step = array_sz / no_buckets, range_end = range_step;
        uint64_t range_mod = array_sz % no_buckets;
        for(uint64_t i = 0; i < no_buckets; i++){
            if(i < range_mod) range_end++;
            tasks.push_back( async(launch::async, create_partition, i, range_start, range_end) );
            range_start = range_end;
            range_end += range_step;
        }
        // wait for all tasks to finish
        for(auto& t: tasks) t.get();
    }

    { // compute [sequentially] the size of each local permutation
        uint64_t offset = 0;
        for (uint64_t bucket_id = 0, sz = buckets.size(); bucket_id < sz; bucket_id++) {
            Bucket<T>* bucket = buckets[bucket_id].get();
            uint64_t size = 0;
            assert((size_t) bucket->m_chunks_size == buckets.size());
            for(int i = 0; i < bucket->m_chunks_size; i++){ size += bucket->m_chunks[i]->size(); }

            bucket->m_permutation = array + offset;
            bucket->m_permutation_sz = size;
            offset += size;
        }
    }

    { // perform a local permutation
        auto local_permutation = [&buckets/*, bytes_per_element*/](int bucket_id){
            // store all the values in a single array
            Bucket<T>* bucket = buckets[bucket_id].get();
            T* __restrict permutation = bucket->m_permutation; //new CByteArray(bytes_per_element, capacity);
            size_t index = 0;
            for(int i = 0; i < bucket->m_chunks_size; i++){
                vector<T>* vector = bucket->m_chunks[i];
                for(size_t j = 0, sz = vector->size(); j < sz; j++){
//                  bucket->permutation->set_value_at(index++, vector->at(j));
                    permutation[index++] = vector->at(j);
                }
            }

            // deallocate the chunks to save some memory
            for(int i = 0; i < bucket->m_chunks_size; i++){
                delete bucket->m_chunks[i]; bucket->m_chunks[i] = nullptr;
            }
            delete[] bucket->m_chunks; bucket->m_chunks = nullptr;

            // perform a local permutation
            if(bucket->m_permutation_sz >= 2){
                T* __restrict permutation = bucket->m_permutation;
                for(size_t i = 0; i < bucket->m_permutation_sz -2; i++){
                    size_t j = uniform_int_distribution<size_t>{i, bucket->m_permutation_sz -1}(bucket->m_random_generator);
                    std::swap(permutation[i], permutation[j]);  // swap A[i] with A[j]
                }
            }
        };

        std::vector<future<void>> tasks;
        for(size_t i = 0; i < no_buckets; i++){
            tasks.push_back( async(launch::async, local_permutation, i) );
        }
        // wait for all tasks to finish
        for(auto& t: tasks) t.get();
    }

    // there is no need to explicitly concatenate all the results together as before (CByteArray::merge)
    // done
}


template<typename T>
void permute(T* array, uint64_t array_sz, uint64_t seed){
    // trampoline to the implementation
    do_permute(array, array_sz, /* no_buckets = */ max(4u, thread::hardware_concurrency()) * 8, seed);
}

// explicit instantiations
template void permute<int32_t>(int32_t* array, uint64_t, uint64_t seed);
template void permute<uint32_t>(uint32_t* array, uint64_t, uint64_t seed);
template void permute<int64_t>(int64_t* array, uint64_t, uint64_t seed);
template void permute<uint64_t>(uint64_t* array, uint64_t, uint64_t seed);

}