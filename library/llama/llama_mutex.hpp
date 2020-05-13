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

#include <future>
#include <pthread.h>
#include <utility>

#include "common/circular_array.hpp"

namespace gfe::library::llama_details {

/**
 * A simple implementation of a fair read/write lock. This is used in LLAMA to guarantee
 * that the thread involved in the checkpoint (compactation) does not suffer starvation.
 */
class FairSharedMutex {
    FairSharedMutex(const FairSharedMutex&) = delete;
    FairSharedMutex& operator=(const FairSharedMutex&) = delete;

    pthread_spinlock_t m_lock; // spin lock
    int m_state; // -1: write, 0: free, > 0 num readers
    struct Item {
        int m_type; // -1: writer, +1: reader
        std::promise<void>* m_producer; // waiting thread
    };
    common::CircularArray<Item> m_waitlist; // threads waiting to acquire the lock

    // Wait the next turn to acquire the mutex
    void wait(int role /* -1 => writer, +1 => reader */);

    // Wake the next item in the queue
    void wake_next();

public:
    // Constructor
    FairSharedMutex();

    // Destructor
    ~FairSharedMutex();

    // Exclusively lock the mutex
    void lock();

    // Try to exclusively lock the mutex
    bool try_lock();

    // Release the exclusive lock
    void unlock();

    // Acquire the mutex as a reader
    void lock_shared();

    // Try to acquire the shared lock and return immediately, that is, don't block.
    bool try_lock_shared();

    // Release the mutex as a reader
    void unlock_shared();
};

/*****************************************************************************
 *                                                                           *
 *   Implementation details                                                  *
 *                                                                           *
 *****************************************************************************/
inline
FairSharedMutex::FairSharedMutex() : m_state(0) {
    [[maybe_unused]] int rc = pthread_spin_init(&m_lock, /* flags */ 0);
    assert(rc == 0 && "Cannot initialise the spin lock");
}

inline
FairSharedMutex::~FairSharedMutex() {
    pthread_spin_destroy(&m_lock);
}

inline
void FairSharedMutex::lock(){
    pthread_spin_lock(&m_lock);
    while(m_state != 0){ wait(/* writer */ -1); }
    m_state = -1;
    pthread_spin_unlock(&m_lock);
}

inline
bool FairSharedMutex::try_lock(){
    bool is_locked = false;
    pthread_spin_lock(&m_lock);
    if(m_state == 0){
        m_state = -1;
        is_locked = true;
    }
    pthread_spin_unlock(&m_lock);
    return is_locked;
}

inline
void FairSharedMutex::unlock(){
    pthread_spin_lock(&m_lock);
    assert(m_state == -1 && "It should have been previously locked");
    m_state = 0;
    wake_next();
    pthread_spin_unlock(&m_lock);
}

inline
void FairSharedMutex::lock_shared(){
    pthread_spin_lock(&m_lock);
    if(m_state == 0 || (m_state > 0 && m_waitlist.empty())){
        m_state ++;
        pthread_spin_unlock(&m_lock);
        return;
    }

    do { wait( /* reader */ +1); } while (m_state < 0);

    m_state++;
    pthread_spin_unlock(&m_lock);
}

inline
bool FairSharedMutex::try_lock_shared(){
    bool is_locked = false;
    pthread_spin_lock(&m_lock);
    if(m_state == 0 || (m_state > 0 && m_waitlist.empty())){
        m_state ++;
        is_locked = true;
    }
    pthread_spin_unlock(&m_lock);

    return is_locked;
}

inline
void FairSharedMutex::unlock_shared() {
    pthread_spin_lock(&m_lock);
    assert(m_state > 0 && "The shared lock should have been previously acquired");
    m_state--;
    if(m_state == 0) wake_next();
    pthread_spin_unlock(&m_lock);
}

inline
void FairSharedMutex::wait(int role){
    assert(role == -1 /* writer */ || role == +1 /* reader */);
    std::promise<void> producer;
    std::future<void> consumer = producer.get_future();
    m_waitlist.append(Item{role, &producer});
    pthread_spin_unlock(&m_lock);
    consumer.get();
    pthread_spin_lock(&m_lock);
}

inline
void FairSharedMutex::wake_next() {
    if(m_waitlist.empty()) return;
    switch(m_waitlist[0].m_type){
    case -1: // writer
        m_waitlist[0].m_producer->set_value();
        m_waitlist.pop();
        break;
    case +1: // reader
        do { // wake the sequence of readers
            m_waitlist[0].m_producer->set_value();
            m_waitlist.pop();
        } while(!m_waitlist.empty() && m_waitlist[0].m_type == +1);
        break;
    default:
        assert(0 && "Invalid type");
    }
}


} // namespace

