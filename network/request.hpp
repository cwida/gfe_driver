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

#include <cinttypes>
#include <ostream>

namespace network {


/**
 * Type of message request
 */
enum class RequestType : uint32_t {
    TERMINATE_SERVER, // terminate the server, imply this connection as well
    TERMINATE_WORKER, // terminate only the worker connection
    ON_MAIN_INIT, ON_THREAD_INIT, ON_THREAD_DESTROY, ON_MAIN_DESTROY,
    NUM_EDGES, NUM_VERTICES,
    ADD_VERTEX, DELETE_VERTEX, ADD_EDGE, DELETE_EDGE
};

/**
 * Generic request
 */
class Request {
public:
    uint32_t m_message_size; // the whole size of the request message, in bytes
    RequestType m_type; // the type of request
};

/**
 * Send a request such as on_main_init, on_thread_init, on_thread_destroy, on_main_destroy
 */
template<int N>
class GenericRequest : public Request {
public:
    uint64_t m_arguments[N];
};

/**
 * Dump to the output stream the content of the request
 */
std::ostream& operator<<(std::ostream& out, const Request& request);
std::ostream& operator<<(std::ostream& out, const Request* request);
std::ostream& operator<<(std::ostream& out, RequestType type);

}


