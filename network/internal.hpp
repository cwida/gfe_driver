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

#include "network/error.hpp"

#include <cerrno>
#include <mutex>
#include <string.h> // strerror

// Handy macro to throw a network error with ERROR( message )
#undef CURRENT_ERROR_TYPE
#define CURRENT_ERROR_TYPE ::gfe::network::NetworkError

// Macro to report an error with the content of errno
#define ERROR_ERRNO(message) ERROR(message << ". Low level description: " << strerror(errno) << " (errno: " << errno << ")")

namespace gfe::network {

// Some functions from the C library are not thread-safe, such as gethostbyname() and gethostbyaddr().
// More details: http://man7.org/linux/man-pages/man3/gethostbyname.3.html
// We wrap each call to these functions inside the following global lock
extern std::mutex g_network_lock;

// Convenience macro to wrap a critical section inside g_network_lock
#define SYNCHRONISE_NETWORK(stmt) { std::scoped_lock<std::mutex> lock(g_network_lock); stmt; }

} // namespace
