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
    HAS_VERTEX, HAS_EDGE,
    LOAD, // load the graph from disk
    ADD_VERTEX, DELETE_VERTEX, ADD_EDGE, DELETE_EDGE
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
    const uint32_t m_message_size;
    const Type m_type;

public:
    template<typename... Args>
    Message(Type type, Args... args);

    Message(Type type, const char* string);

    // Retrieve the total number of arguments in the message
    int num_arguments() const;

    // Retrieve the size of the message
    size_t message_size() const;

    // The type associated to this message
    Type type() const;

    // Get the given argument
    template<typename T = uint64_t>
    T get(int index) const;

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
static typename std::enable_if_t< std::is_signed_v<TSigned> >
store_single_arg(char* buffer, TSigned arg){
    reinterpret_cast<int64_t*>(buffer)[0] = arg;
}

template<typename TUnsigned>
static typename std::enable_if_t< std::is_unsigned_v<TUnsigned> >
store_single_arg(char* buffer, TUnsigned arg){
    reinterpret_cast<uint64_t*>(buffer)[0] = arg;
}

static inline void store_args(char* buffer){
    // no args to store...
}

// base case of the recursion
template<typename T>
static void store_args(char* buffer, T arg){
    store_single_arg(buffer, arg);
}

template<typename T, typename... Args>
static void store_args(char* buffer, T first, Args... rest){
    store_single_arg(buffer, first);
    store_args(buffer + sizeof(uint64_t), rest...);
}

template<typename TSigned>
static typename std::enable_if_t< std::is_signed_v<TSigned>, TSigned > retrieve_single_arg(const char* buffer){
    return static_cast<TSigned>(reinterpret_cast<const int64_t*>(buffer)[0]);
}

template<typename TUnsigned>
static typename std::enable_if_t< std::is_unsigned_v<TUnsigned>, TUnsigned > retrieve_single_arg(const char* buffer){
    return static_cast<TUnsigned>(reinterpret_cast<const uint64_t*>(buffer)[0]);
}

} // namespace details

// Implementation details
template<typename Type>
template<typename... Args>
Message<Type>::Message(Type type, Args... args) : m_message_size(sizeof(Message<Type>) + sizeof...(args) * sizeof(uint64_t)), m_type(type){
    details::store_args(buffer(), args...);
}

template<typename Type>
Message<Type>::Message(Type type, const char* string) : m_message_size(sizeof(Message<Type>) + strlen(string) +1), m_type(type){
    strcpy(buffer(), string); // a bit of hack
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



