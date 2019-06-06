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

#include "response.hpp"

using namespace std;

namespace network {

Response::Response(uint32_t message_size, ResponseType type) : m_message_size(message_size), m_type(type) { }
ResponseSingleValue::ResponseSingleValue(uint32_t message_size, ResponseType type, uint64_t return_value) :
        Response(message_size, type), m_return_value(return_value) { }

std::ostream& operator<<(std::ostream& out, ResponseType type){
    switch(type){
    case ResponseType::OK: out << "OK"; break;
    case ResponseType::SINGLE_VALUE: out << "SINGLE_VALUE"; break;
    default: out << "UNKNOWN (response code: " << (uint32_t) type << ")";
    }
    return out;
}

std::ostream& operator<<(std::ostream& out, const Response& response){
    out << "[RESPONSE " << response.m_type << ", message size: " << response.m_message_size;

    switch(response.m_type){
    case ResponseType::SINGLE_VALUE: {
        out << ", value: " << reinterpret_cast<const ResponseSingleValue*>(&response)->m_return_value;
    } break;
    default:
        ; /* nop */
    }

    out << "]";
    return out;
}

std::ostream& operator<<(std::ostream& out, const Response* response){
    if(response == nullptr){
        out << "[RESPONSE nullptr]";
    } else {
        out << *response;
    }
    return out;
}

} // namespace network
