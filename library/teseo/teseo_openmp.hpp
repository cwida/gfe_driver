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

#include <sched.h> // cpu_set_t
#include "teseo.hpp"

namespace gfe::library { class TeseoDriver; } // forward decl.

namespace gfe::library::teseo_driver_internal {

/**
 * Restrict the thread affinity to one of the sockets in the machine
 */
class ThreadAffinity {
    ThreadAffinity(const ThreadAffinity&) = delete; // thread specific
    ThreadAffinity& operator=(const ThreadAffinity&) = delete;

    const bool m_is_enabled;
    cpu_set_t* m_copy_to_restore = nullptr; // != nullptr => on dtor reset the thread affinity of the thread

public:
    /**
     * Restrict the thread affinity to one of the sockets of the machine
     */
    ThreadAffinity(bool enabled);

    /**
     * Restore the thread affinity of the thread
     */
    ~ThreadAffinity();

    /**
     * Check whether the thread affinity has been set in the current thread
     */
    bool is_enabled() const noexcept;
};

/**
 * RAII to register/unregister by an OpenMP worker thread to Teseo. A ID of a worker thread
 * in OpenMP is always greater than 0.
 *
 * All threads used in Teseo must be registered for the correct execution of its Garbage Collector.
 */
class RegisterThread {
    RegisterThread& operator=(const RegisterThread&) = delete;
    uint64_t m_encoded_addr; // The address of Teseo + 1 bit to check whether to unregister the thread

    // Teseo address
    ::teseo::Teseo* teseo() const;

    // Check whether to unregister the thread
    bool is_enabled() const;

public:
    /**
     * Register the thread to Teseo
     */
    RegisterThread(::teseo::Teseo* teseo);

    /**
     * Copy constructor
     */
    RegisterThread(const RegisterThread& rt);

    /**
     * Destructor, unregister the thread
     */
    ~RegisterThread();
};


/**
 * All machinery needed to properly register/unregister the OpenMP threads in Teseo
 */
class OpenMP {
    OpenMP& operator=(const OpenMP&) = delete;

    ThreadAffinity m_thread_affinity; // restrict the execution of the threads in a given socket
    RegisterThread m_thread_tracking; // register/unregister the thread in Teseo
    teseo::Transaction m_transaction; // transaction in use
    teseo::Iterator m_iterator; // iterator in use

public:
    /**
     * Initialise the OpenMP state for the master thread
     */
    OpenMP(TeseoDriver* driver);

    /**
     * Initialise the OpenMP state for the worker threads
     */
    OpenMP(const OpenMP& openmp);

    /**
     * Destructor
     */
    ~OpenMP();

    /**
     * Retrieve the underlying transaction
     */
    teseo::Transaction& transaction() noexcept {
        return m_transaction;
    }

    /**
     * Retrieve the underlying iterator
     */
    teseo::Iterator& iterator() noexcept {
        return m_iterator;
    }
};


} // namespace
