// The following source code has been extracted from the GAP Benchmark Suite
// https://github.com/sbeamer/gapbs
// The author is Scott Beamer
//
// Copyright (c) 2015, The Regents of the University of California (Regents)
// All Rights Reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
// 3. Neither the name of the Regents nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL REGENTS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#pragma once

#include <algorithm>
#include <cinttypes>
#include <cstddef> // size_t

namespace gapbs {

/*
GAP Benchmark Suite
File:   Platform Atomics
Author: Scott Beamer
Wrappers for compiler intrinsics for atomic memory operations (AMOs)
 - If not using OpenMP (serial), provides serial fallbacks
*/

#if defined _OPENMP

#if defined __GNUC__

// gcc/clang/icc instrinsics

template<typename T, typename U>
T fetch_and_add(T &x, U inc) {
    return __sync_fetch_and_add(&x, inc);
}

template<typename T>
bool compare_and_swap(T &x, const T &old_val, const T &new_val) {
    return __sync_bool_compare_and_swap(&x, old_val, new_val);
}

template<>
inline bool compare_and_swap(float &x, const float &old_val, const float &new_val) {
    return __sync_bool_compare_and_swap(reinterpret_cast<uint32_t*>(&x),
            reinterpret_cast<const uint32_t&>(old_val),
            reinterpret_cast<const uint32_t&>(new_val));
}

template<>
inline bool compare_and_swap(double &x, const double &old_val, const double &new_val) {
    return __sync_bool_compare_and_swap(reinterpret_cast<uint64_t*>(&x),
            reinterpret_cast<const uint64_t&>(old_val),
            reinterpret_cast<const uint64_t&>(new_val));
}

#elif __SUNPRO_CC

// sunCC (solaris sun studio) intrinsics
// less general, only work for int32_t, int64_t, uint32_t, uint64_t
// http://docs.oracle.com/cd/E19253-01/816-5168/6mbb3hr06/index.html

#include <atomic.h>
#include <cinttypes>

int32_t fetch_and_add(int32_t &x, int32_t inc) {
    return atomic_add_32_nv((volatile uint32_t*) &x, inc) - inc;
}

int64_t fetch_and_add(int64_t &x, int64_t inc) {
    return atomic_add_64_nv((volatile uint64_t*) &x, inc) - inc;
}

uint32_t fetch_and_add(uint32_t &x, uint32_t inc) {
    return atomic_add_32_nv((volatile uint32_t*) &x, inc) - inc;
}

uint64_t fetch_and_add(uint64_t &x, uint64_t inc) {
    return atomic_add_64_nv((volatile uint64_t*) &x, inc) - inc;
}

bool compare_and_swap(int32_t &x, const int32_t &old_val, const int32_t &new_val) {
    return old_val == atomic_cas_32((volatile uint32_t*) &x, old_val, new_val);
}

bool compare_and_swap(int64_t &x, const int64_t &old_val, const int64_t &new_val) {
    return old_val == atomic_cas_64((volatile uint64_t*) &x, old_val, new_val);
}

bool compare_and_swap(uint32_t &x, const uint32_t &old_val, const uint32_t &new_val) {
    return old_val == atomic_cas_32((volatile uint32_t*) &x, old_val, new_val);
}

bool compare_and_swap(uint64_t &x, const uint64_t &old_val, const uint64_t &new_val) {
    return old_val == atomic_cas_64((volatile uint64_t*) &x, old_val, new_val);
}

bool compare_and_swap(float &x, const float &old_val, const float &new_val) {
    return old_val == atomic_cas_32((volatile uint32_t*) &x,
            (const volatile uint32_t&) old_val,
            (const volatile uint32_t&) new_val);
}

bool compare_and_swap(double &x, const double &old_val, const double &new_val) {
    return old_val == atomic_cas_64((volatile uint64_t*) &x,
            (const volatile uint64_t&) old_val,
            (const volatile uint64_t&) new_val);
}

#else   // defined __GNUC__ __SUNPRO_CC

#error No atomics available for this compiler but using OpenMP

#endif  // else defined __GNUC__ __SUNPRO_CC

#else   // defined _OPENMP

// serial fallbacks

template<typename T, typename U>
T fetch_and_add(T &x, U inc) {
    T orig_val = x;
    x += inc;
    return orig_val;
}

template<typename T>
bool compare_and_swap(T &x, const T &old_val, const T &new_val) {
    if (x == old_val) {
        x = new_val;
        return true;
    }
    return false;
}

#endif  // else defined _OPENMP

/*
GAP Benchmark Suite
Class:  Bitmap
Author: Scott Beamer
Parallel bitmap that is thread-safe
 - Can set bits in parallel (set_bit_atomic) unlike std::vector<bool>
 */

class Bitmap {
public:
    explicit Bitmap(size_t size) {
        uint64_t num_words = (size + kBitsPerWord - 1) / kBitsPerWord;
        start_ = new uint64_t[num_words];
        end_ = start_ + num_words;
    }

    ~Bitmap() {
        delete[] start_;
    }

    void reset() {
        std::fill(start_, end_, 0);
    }

    void set_bit(size_t pos) {
        start_[word_offset(pos)] |= ((uint64_t) 1l << bit_offset(pos));
    }

    void set_bit_atomic(size_t pos) {
        uint64_t old_val, new_val;
        do {
            old_val = start_[word_offset(pos)];
            new_val = old_val | ((uint64_t) 1l << bit_offset(pos));
        } while (!compare_and_swap(start_[word_offset(pos)], old_val, new_val));
    }

    bool get_bit(size_t pos) const {
        return (start_[word_offset(pos)] >> bit_offset(pos)) & 1l;
    }

    void swap(Bitmap &other) {
        std::swap(start_, other.start_);
        std::swap(end_, other.end_);
    }

private:
    uint64_t *start_;
    uint64_t *end_;

    static const uint64_t kBitsPerWord = 64;
    static uint64_t word_offset(size_t n) { return n / kBitsPerWord; }
    static uint64_t bit_offset(size_t n) { return n & (kBitsPerWord - 1); }
};

/*
GAP Benchmark Suite
Class:  SlidingQueue
Author: Scott Beamer
Double-buffered queue so appends aren't seen until SlideWindow() called
 - Use QueueBuffer when used in parallel to avoid false sharing by doing
   bulk appends from thread-local storage
*/

template <typename T>
class QueueBuffer;

template <typename T>
class SlidingQueue {
    T *shared;
    size_t shared_in;
    size_t shared_out_start;
    size_t shared_out_end;
    friend class QueueBuffer<T>;

public:
    explicit SlidingQueue(size_t shared_size) {
        shared = new T[shared_size];
        reset();
    }

    ~SlidingQueue() {
        delete[] shared;
    }

    void push_back(T to_add) {
        shared[shared_in++] = to_add;
    }

    bool empty() const {
        return shared_out_start == shared_out_end;
    }

    void reset() {
        shared_out_start = 0;
        shared_out_end = 0;
        shared_in = 0;
    }

    void slide_window() {
        shared_out_start = shared_out_end;
        shared_out_end = shared_in;
    }

    typedef T* iterator;

    iterator begin() const {
        return shared + shared_out_start;
    }

    iterator end() const {
        return shared + shared_out_end;
    }

    size_t size() const {
        return end() - begin();
    }
};

template <typename T>
class QueueBuffer {
    size_t in;
    T *local_queue;
    SlidingQueue<T> &sq;
    const size_t local_size;

public:
    explicit QueueBuffer(SlidingQueue<T> &master, size_t given_size = 16384)
    : sq(master), local_size(given_size) {
        in = 0;
        local_queue = new T[local_size];
    }

    ~QueueBuffer() {
        delete[] local_queue;
    }

    void push_back(T to_add) {
        if (in == local_size)
            flush();
        local_queue[in++] = to_add;
    }

    void flush() {
        T *shared_queue = sq.shared;
        size_t copy_start = fetch_and_add(sq.shared_in, in);
        std::copy(local_queue, local_queue+in, shared_queue+copy_start);
        in = 0;
    }
};

/*
GAP Benchmark Suite
Class:  pvector
Author: Scott Beamer

Vector class with ability to not initialize or do initialize in parallel
 - std::vector (when resizing) will always initialize, and does it serially
 - When pvector is resized, new elements are uninitialized
 - Resizing is not thread-safe
*/

template <typename T_>
class pvector {
public:
    typedef T_* iterator;

    pvector() : start_(nullptr), end_size_(nullptr), end_capacity_(nullptr) {}

    explicit pvector(size_t num_elements) {
        start_ = new T_[num_elements];
        end_size_ = start_ + num_elements;
        end_capacity_ = end_size_;
    }

    pvector(size_t num_elements, T_ init_val) : pvector(num_elements) {
        fill(init_val);
    }

    pvector(iterator copy_begin, iterator copy_end)
    : pvector(copy_end - copy_begin) {
#pragma omp parallel for
        for (size_t i=0; i < capacity(); i++)
            start_[i] = copy_begin[i];
    }

    // don't want this to be copied, too much data to move
    pvector(const pvector &other) = delete;

    // prefer move because too much data to copy
    pvector(pvector &&other)
    : start_(other.start_), end_size_(other.end_size_),
      end_capacity_(other.end_capacity_) {
        other.start_ = nullptr;
        other.end_size_ = nullptr;
        other.end_capacity_ = nullptr;
    }

    // want move assignment
    pvector& operator= (pvector &&other) {
        start_ = other.start_;
        end_size_ = other.end_size_;
        end_capacity_ = other.end_capacity_;
        other.start_ = nullptr;
        other.end_size_ = nullptr;
        other.end_capacity_ = nullptr;
        return *this;
    }

    ~pvector() {
        if (start_ != nullptr)
            delete[] start_;
    }

    // not thread-safe
    void reserve(size_t num_elements) {
        if (num_elements > capacity()) {
            T_ *new_range = new T_[num_elements];
#pragma omp parallel for
            for (size_t i=0; i < size(); i++)
                new_range[i] = start_[i];
            end_size_ = new_range + size();
            delete[] start_;
            start_ = new_range;
            end_capacity_ = start_ + num_elements;
        }
    }

    bool empty() {
        return end_size_ == start_;
    }

    void clear() {
        end_size_ = start_;
    }

    void resize(size_t num_elements) {
        reserve(num_elements);
        end_size_ = start_ + num_elements;
    }

    T_& operator[](size_t n) {
        return start_[n];
    }

    const T_& operator[](size_t n) const {
        return start_[n];
    }

    void push_back(T_ val) {
        if (size() == capacity()) {
            size_t new_size = capacity() == 0 ? 1 : capacity() * growth_factor;
            reserve(new_size);
        }
        *end_size_ = val;
        end_size_++;
    }

    void fill(T_ init_val) {
#pragma omp parallel for
        for (T_* ptr=start_; ptr < end_size_; ptr++)
            *ptr = init_val;
    }

    size_t capacity() const {
        return end_capacity_ - start_;
    }

    size_t size() const {
        return end_size_ - start_;
    }

    iterator begin() const {
        return start_;
    }

    iterator end() const {
        return end_size_;
    }

    T_* data() const {
        return start_;
    }

    void swap(pvector &other) {
        std::swap(start_, other.start_);
        std::swap(end_size_, other.end_size_);
        std::swap(end_capacity_, other.end_capacity_);
    }


private:
    T_* start_;
    T_* end_size_;
    T_* end_capacity_;
    static const size_t growth_factor = 2;
};

} // namespace
