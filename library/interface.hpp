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

#include <iostream>

namespace library {

/**
 * Base interface, implemented by all systems
 */
class Interface {
public:
    /**
     * Dummy constructor
     */
    Interface();

    /**
     * Virtual destructor
     */
    virtual ~Interface();

    /**
     * Dump the content of the graph to stdout
     */
    virtual void dump() const = 0;

    /**
     * Thread initialisation callbacks
     */
    // Invoked at the start of the experiment by the controller thread with the number of threads that will be used
    virtual void on_main_init(int num_threads);
    // Invoked by each worker thread separately, with a unique thread id, in [0, num_threads)
    virtual void on_thread_init(int thread_id);
    // Invoked by each worker thread separately, with the same thread id given at on_thread_init
    virtual void on_thread_destroy(int thread_id);
    // Invoked at the end of the experiment by the controller thread
    virtual void on_main_destroy();
};

} // namespace library

