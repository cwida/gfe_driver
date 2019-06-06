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
 * The type of the message sent back to the client
 */
enum class ResponseType: uint32_t {
    OK, // it doesn't contain any return value, usually done for RPC that do not have a return value
    SINGLE_VALUE, // it contains a single return value of 8 bytes
};

/**
 * The response message sent from the server to the client
 */
class Response {
public:
    uint32_t m_message_size; // the whole size of the request message, in bytes
    ResponseType m_type; // the type of request

    Response(uint32_t message_size, ResponseType type);
};

/**
 * A response with a single return value
 */
class ResponseSingleValue : public Response {
public:
    uint64_t m_return_value; // a single return value of 8 bytes (int64_t, uint64_t)

    ResponseSingleValue(uint32_t message_size, ResponseType type, uint64_t return_value);
};

/**
 * Dump to the output stream the content of the request
 */
std::ostream& operator<<(std::ostream& out, const Response& response);
std::ostream& operator<<(std::ostream& out, const Response* response);
std::ostream& operator<<(std::ostream& out, ResponseType type);


} // namespace network
