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

#include <cassert>
#include <cinttypes>
#include <cmath>
#include <cstring>
#include <ostream>
#include <type_traits>

namespace network {

/**
 * Type of request message, sent from the client to the server
 */
enum class RequestType : uint32_t {
    TERMINATE_SERVER, // terminate the server, imply this connection as well
    TERMINATE_WORKER, // terminate only the worker connection
    ON_MAIN_INIT, ON_THREAD_INIT, ON_THREAD_DESTROY, ON_MAIN_DESTROY,
    NUM_EDGES, NUM_VERTICES,
    HAS_VERTEX, HAS_EDGE, GET_WEIGHT,
    LOAD, // load the graph from disk
    ADD_VERTEX, DELETE_VERTEX, ADD_EDGE, DELETE_EDGE,
    DUMP_STDOUT, DUMP_FILE, // #dump()
    BFS, PAGERANK, WCC, CDLP, LCC, SSSP // graphalytics interface
};

/**
 * The type of a response message, sent back from the server to the client
 */
enum class ResponseType: uint32_t {
    OK, // ok, the request has been processed correctly
    NOT_SUPPORTED, // the remote server does not support the given operation
};

/**
 * A generic message, type + arguments, sent between the clients and the server
 */
template<typename Type>
class Message {
    uint32_t m_message_size;
    const Type m_type;

public:
    template<typename... Args>
    Message(Type type, Args... args);

    // Retrieve the total number of arguments in the message
    int num_arguments() const;

    // Retrieve the size of the message
    size_t message_size() const;

    // The type associated to this message
    Type type() const;

    // Get the given argument
    template<typename T = uint64_t>
    T get(int index) const;

    // Get the given argument as a string
    std::string get_string(int index) const;

    // Get the space where the arguments are stored
    const char* buffer() const;
    char* buffer();
};


// Specific message sent from the clients to the server
using Request = Message<RequestType>;

// Specific message sent from the server to the client
using Response = Message<ResponseType>;

// Store multiple arguments in the given buffer
namespace details {
template<typename TSigned>
static typename std::enable_if_t< std::is_integral_v<TSigned> &&std::is_signed_v<TSigned>, uint64_t>
store_single_arg(char* buffer, TSigned arg){
    reinterpret_cast<int64_t*>(buffer)[0] = arg;
    return sizeof(int64_t); // size of the argument
}

template<typename TUnsigned>
static typename std::enable_if_t< std::is_integral_v<TUnsigned> && std::is_unsigned_v<TUnsigned>, uint64_t >
store_single_arg(char* buffer, TUnsigned arg){
    reinterpret_cast<uint64_t*>(buffer)[0] = arg;
    return sizeof(uint64_t); // size of the argument
}

template<typename TFloat>
static typename std::enable_if_t< std::is_floating_point_v<TFloat> && sizeof(TFloat) <= sizeof(double), uint64_t >
store_single_arg(char* buffer, TFloat arg){
    reinterpret_cast<double*>(buffer)[0] = arg;
    return sizeof(double); // size of the argument
}

inline static uint64_t
store_single_arg(char* buffer, const char* str){
    uint64_t* length = reinterpret_cast<uint64_t*>(buffer);
    if(str == nullptr){
        *length = 0;
    } else {
        *length = strlen(str);
        strcpy(buffer + sizeof(uint64_t), str); // can overflow
    }

    return sizeof(uint64_t) + ceil(static_cast<double>(*length) / sizeof(uint64_t)); // round up
}

inline static uint64_t
store_single_arg(char* buffer, const std::string& str){
    return store_single_arg(buffer, str.c_str());
}

static inline uint64_t store_args(char* buffer){
    return 0; // no args to store...
}

// base case of the recursion
template<typename T>
static uint64_t store_args(char* buffer, T arg){
    return store_single_arg(buffer, arg);
}

template<typename T, typename... Args>
static uint64_t store_args(char* buffer, T first, Args... rest){
    uint64_t offset = store_single_arg(buffer, first);
    return offset + store_args(buffer + offset, rest...);
}

template<typename TSigned>
static typename std::enable_if_t< std::is_integral_v<TSigned> && std::is_signed_v<TSigned>, TSigned > retrieve_single_arg(const char* buffer){
    return static_cast<TSigned>(reinterpret_cast<const int64_t*>(buffer)[0]);
}

template<typename TUnsigned>
static typename std::enable_if_t< std::is_integral_v<TUnsigned> && std::is_unsigned_v<TUnsigned>, TUnsigned > retrieve_single_arg(const char* buffer){
    return static_cast<TUnsigned>(reinterpret_cast<const uint64_t*>(buffer)[0]);
}

template<typename TFloat>
static typename std::enable_if_t< std::is_floating_point_v<TFloat>, TFloat > retrieve_single_arg(const char* buffer){
    return static_cast<TFloat>(reinterpret_cast<const double*>(buffer)[0]);
}

} // namespace details

// Implementation details
template<typename Type>
template<typename... Args>
Message<Type>::Message(Type type, Args... args) : m_message_size(sizeof(Message<Type>)), m_type(type){
    m_message_size += details::store_args(buffer(), args...);
}


template<typename Type>
int Message<Type>::num_arguments() const {
    return (static_cast<int>(message_size()) - sizeof(Message<Type>)) / sizeof(uint64_t);
}

template<typename Type>
size_t Message<Type>::message_size() const {
    return m_message_size;
}

template<typename Type>
Type Message<Type>::type() const {
    return m_type;
}

template<typename Type>
const char* Message<Type>::buffer() const{
    return reinterpret_cast<const char*>(this) + sizeof(Message<Type>);
}

template<typename Type>
char* Message<Type>::buffer(){
    return reinterpret_cast<char*>(this) + sizeof(Message<Type>);
}

template<typename Type>
template<typename T>
T Message<Type>::get(int index) const{
    assert(index >= 0 && "Invalid index, negative value");
    assert(index < num_arguments() && "Out of bound");
    return details::retrieve_single_arg<T>(buffer() + sizeof(uint64_t) * index);
}

template<typename Type>
std::string Message<Type>::get_string(int index) const {
    assert(index >= 0 && "Invalid index, negative value");
    assert(index < num_arguments() && "Out of bound");
    auto length = details::retrieve_single_arg<uint64_t>(buffer() + sizeof(uint64_t) * index);
    return std::string { buffer() + sizeof(uint64_t) * (index+1), length };
}

/**
 * Dump to the output stream the content of a message
 */
std::ostream& operator<<(std::ostream& out, const Request& request);
std::ostream& operator<<(std::ostream& out, const Request* request);
std::ostream& operator<<(std::ostream& out, RequestType type);
std::ostream& operator<<(std::ostream& out, const Response& response);
std::ostream& operator<<(std::ostream& out, const Response* response);
std::ostream& operator<<(std::ostream& out, ResponseType type);

} // namespace network



