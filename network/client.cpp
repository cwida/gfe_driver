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

#include "client.hpp"

#include <arpa/inet.h>
#include <cassert>
#include <cstring>
#include <netdb.h> // gethostbyname
#include <sys/socket.h>
#include <type_traits>
#include <unistd.h>

#include "configuration.hpp"
#include "internal.hpp"

using namespace std;

namespace network {

#define RPC_ERROR(message) RAISE_EXCEPTION(RPCError, "Error propagated by the remote server:\n" << message)
#define TIMEOUT_ERROR RAISE_EXCEPTION(library::TimeoutError, "The operation requested timed out");

/*****************************************************************************
 *                                                                           *
 * Initialisation                                                            *
 *                                                                           *
 *****************************************************************************/

thread_local int Client::m_worker_id { 0 };

Client::Client(const std::string& host, int port) : m_server_host(host), m_server_port(port) {
    // reset the content of the connections
    for(int i = 0; i < max_num_connections; i++){
        m_connections[i].m_fd = -1;
        m_connections[i].m_buffer_read = m_connections[i].m_buffer_write = nullptr;
    }

//    LOG("[client] Connecting to " << m_server_host << ":" << m_server_port << " ...");
    connect();
    LOG("[client] Connection established!");
}

Client::~Client(){
    m_worker_id = 0;
    request(RequestType::TERMINATE_WORKER);

    // close all connections still open
    for(int i = 0; i < max_num_connections; i++){ disconnect(i); }
}

void Client::connect(){
    if(m_connections[m_worker_id].m_fd >= 0) return; // already connected

    int fd = socket(AF_INET, SOCK_STREAM, 0);
    if(fd < 0) ERROR_ERRNO("Cannot initialise the socket");

    // retrieve the address of the host
    struct hostent* remote_address = gethostbyname(m_server_host.c_str());
    if(remote_address == nullptr){
        ERROR("Cannot resolve the address of the remote server `" << m_server_host << "': " << gai_strerror(h_errno));
    }

    struct sockaddr_in address;
    memset(&address, 0, sizeof(address));
    address.sin_family = AF_INET;
    memcpy(&address.sin_addr, remote_address->h_addr, remote_address->h_length);
    address.sin_port = htons(m_server_port);

    int rc = ::connect(fd, (struct sockaddr*) &address, sizeof(address));
    if(rc != 0){
        ERROR_ERRNO("Cannot connect to the remote server " << m_server_host << ":" << m_server_port);
    }

    constexpr uint32_t buffer_default_sz = 4096; // bytes

    m_connections[m_worker_id].m_fd = fd;
    m_connections[m_worker_id].m_buffer_read_sz = buffer_default_sz;
    m_connections[m_worker_id].m_buffer_write_sz = buffer_default_sz;
    m_connections[m_worker_id].m_buffer_read = (char*) malloc(sizeof(char) * buffer_default_sz);
    m_connections[m_worker_id].m_buffer_write = (char*) malloc(sizeof(char) * buffer_default_sz);
}

void Client::disconnect(){
    disconnect(m_worker_id);
}

void Client::disconnect(int worker_id){
    if(worker_id >= max_num_connections) ERROR("Invalid worker_id: " << worker_id);
    if(m_connections[worker_id].m_fd == -1) return;
    close(m_connections[worker_id].m_fd); m_connections[worker_id].m_fd = -1;
    free(m_connections[worker_id].m_buffer_read); m_connections[worker_id].m_buffer_read = nullptr;
    free(m_connections[worker_id].m_buffer_write); m_connections[worker_id].m_buffer_write = nullptr;
}

void Client::terminate_server_on_exit(){
    request(RequestType::TERMINATE_ON_LAST_CONNECTION);
    assert(response()->type() == ResponseType::OK);
}

/*****************************************************************************
 *                                                                           *
 * Request/response                                                          *
 *                                                                           *
 *****************************************************************************/

template<typename... Args>
void Client::request(RequestType type, Args... args){
    assert(m_worker_id >= 0 && m_worker_id < max_num_connections && "Invalid worker id");
    char* buffer = m_connections[m_worker_id].m_buffer_write;
    new (buffer) Request(type, forward<Args>(args)...);
    uint32_t message_sz = reinterpret_cast<uint32_t*>(buffer)[0];
//    cout << "send message_sz: " << message_sz << endl;

    // send the request to the server
    ssize_t bytes_sent = send(m_connections[m_worker_id].m_fd, buffer, message_sz, /* flags */ 0);
    if(bytes_sent == -1) ERROR_ERRNO("send_request, connection error");
    assert(bytes_sent == message_sz && "Message not fully sent");

    // receive the reply from the server
    wait_response();
}

void Client::wait_response() {
    int64_t num_bytes_read { 0 }, recv_bytes { 0 };
    char* buffer = m_connections[m_worker_id].m_buffer_read;
    recv_bytes = recv(m_connections[m_worker_id].m_fd, buffer + num_bytes_read, sizeof(uint32_t), /* flags */ 0);
    if(recv_bytes == -1) ERROR_ERRNO("recv, connection interrupted?");
    assert(recv_bytes == sizeof(uint32_t) && "What the heck have we read?");
    num_bytes_read += recv_bytes;
    uint32_t message_sz = *(reinterpret_cast<uint32_t*>(buffer));
    if(message_sz > m_connections[m_worker_id].m_buffer_read_sz){ // realloc the buffer if it's not large enough
        LOG("realloc buffer, from " << m_connections[m_worker_id].m_buffer_read_sz << " to " << message_sz);
        free(buffer);
        buffer = (char*) malloc(message_sz);
        m_connections[m_worker_id].m_buffer_read = buffer;
        m_connections[m_worker_id].m_buffer_read_sz = message_sz;
        (reinterpret_cast<uint32_t*>(buffer))[0] = message_sz;
    }
    // read the rest of the message
    while(num_bytes_read < message_sz){
        recv_bytes = read(m_connections[m_worker_id].m_fd, buffer + num_bytes_read, message_sz - num_bytes_read);
        if(recv_bytes == -1) ERROR_ERRNO("recv, connection interrupted?");
        num_bytes_read += recv_bytes;
    }
    assert(num_bytes_read == message_sz && "Message read");
}

const Response* Client::response() const {
    assert(m_worker_id >= 0 && m_worker_id < max_num_connections && "Invalid worker id");
    return reinterpret_cast<const Response*>(m_connections[m_worker_id].m_buffer_read);
}

/*****************************************************************************
 *                                                                           *
 * Single RPC calls                                                          *
 *                                                                           *
 *****************************************************************************/

void Client::on_main_init(int num_threads) {
    request(RequestType::ON_MAIN_INIT, num_threads);
    assert(response()->type() == ResponseType::OK);
}

void Client::on_thread_init(int thread_id) {
    if(thread_id < 0 || thread_id > max_num_connections) ERROR("Invalid thread id: " << thread_id);
    m_worker_id = thread_id;
    connect();
    request(RequestType::ON_THREAD_INIT, thread_id);
    assert(response()->type() == ResponseType::OK);
}

void Client::on_thread_destroy(int thread_id) {
    request(RequestType::ON_THREAD_DESTROY, thread_id);
    assert(response()->type() == ResponseType::OK);
    if(thread_id == m_worker_id && m_worker_id > 0){
        request(RequestType::TERMINATE_WORKER);
        disconnect();
    }
    m_worker_id = 0;
}

void Client::on_main_destroy() {
    request(RequestType::ON_MAIN_DESTROY);
    assert(response()->type() == ResponseType::OK);
}

string Client::get_library_name() const {
    const_cast<Client*>(this)->request(RequestType::LIBRARY_NAME);
    assert(response()->type() == ResponseType::OK);
    return response()->get_string(0);
}

uint64_t Client::num_edges() const {
    const_cast<Client*>(this)->request(RequestType::NUM_EDGES);
    assert(response()->type() == ResponseType::OK);
    return response()->get(0);
}

uint64_t Client::num_vertices() const {
    const_cast<Client*>(this)->request(RequestType::NUM_VERTICES);
    assert(response()->type() == ResponseType::OK);
    return response()->get(0);
}

bool Client::is_directed() const {
    const_cast<Client*>(this)->request(RequestType::IS_DIRECTED);
    assert(response()->type() == ResponseType::OK);
    return (bool) response()->get(0);
}

bool Client::has_vertex(uint64_t vertex_id) const {
    const_cast<Client*>(this)->request(RequestType::HAS_VERTEX, vertex_id);
    assert(response()->type() == ResponseType::OK);
    return response()->get<bool>(0);
}

bool Client::has_edge(uint64_t source, uint64_t destination) const {
    const_cast<Client*>(this)->request(RequestType::HAS_EDGE, source, destination);
    assert(response()->type() == ResponseType::OK);
    return response()->get<bool>(0);
}

double Client::get_weight(uint64_t source, uint64_t destination) const {
    const_cast<Client*>(this)->request(RequestType::GET_WEIGHT, source, destination);
    assert(response()->type() == ResponseType::OK);
    return response()->get<double>(0);
}

void Client::load(const std::string& path) {
    request(RequestType::LOAD, path.c_str());
    if(response()->type() == ResponseType::NOT_SUPPORTED){
        ERROR("load(\"" << path << "\"): operation not supported by the remote interface");
    } else if (response()->type() == ResponseType::ERROR){
        RPC_ERROR(response()->get_string(0));
    }
    assert(response()->type() == ResponseType::OK);
}

bool Client::add_vertex(uint64_t vertex_id){
    request(RequestType::ADD_VERTEX, vertex_id);
    if(response()->type() == ResponseType::NOT_SUPPORTED){
        ERROR("add_vertex(" << vertex_id << "): operation not supported by the remote interface");
    }
    assert(response()->type() == ResponseType::OK);
    return response()->get<bool>(0);
}

bool Client::remove_vertex(uint64_t vertex_id){
    request(RequestType::REMOVE_VERTEX, vertex_id);
    if(response()->type() == ResponseType::NOT_SUPPORTED){
        ERROR("delete_vertex(" << vertex_id << "): operation not supported by the remote interface");
    }
    assert(response()->type() == ResponseType::OK);
    return response()->get<bool>(0);
}

bool Client::add_edge(graph::WeightedEdge e){
    request(RequestType::ADD_EDGE, e.source(), e.destination(), e.weight());
    if(response()->type() == ResponseType::NOT_SUPPORTED){
        ERROR("add_edge(" << e.source() << ", " << e.destination() << ", " << e.weight() << "): operation not supported by the remote interface");
    }
    assert(response()->type() == ResponseType::OK);
    return response()->get<bool>(0);
}

bool Client::remove_edge(graph::Edge e){
    request(RequestType::REMOVE_EDGE, e.source(), e.destination());
    if(response()->type() == ResponseType::NOT_SUPPORTED){
        ERROR("delete_edge(" << e.source() << ", " << e.destination() << "): operation not supported by the remote interface");
    }
    assert(response()->type() == ResponseType::OK);
    return response()->get<bool>(0);
}

bool Client::batch(const library::UpdateInterface::SingleUpdate* batch, uint64_t batch_sz, bool force){
    if(batch_sz == 0) return true;
    assert(m_worker_id >= 0 && m_worker_id < max_num_connections && "Invalid worker id");

    constexpr uint32_t header_sz = (uint32_t) sizeof(Request);
    static_assert(header_sz == 8 && "Expected two int32_t fields: the message size and the message type");
    uint32_t body_sz = (uint32_t) sizeof(library::UpdateInterface::SingleUpdate) * batch_sz;
    uint32_t message_sz = header_sz + body_sz;

    if(message_sz > m_connections[m_worker_id].m_buffer_write_sz){
        uint32_t new_size = pow(2, ceil(log2(message_sz))); // next power of 2
        LOG("[worker: " << m_worker_id << "] Reallocate the write buffer to " << new_size << " bytes");
        free(m_connections[m_worker_id].m_buffer_write);
        m_connections[m_worker_id].m_buffer_write = (char*) malloc(new_size);
        m_connections[m_worker_id].m_buffer_write_sz = new_size;
    }

    uint32_t* __restrict buffer_write = reinterpret_cast<uint32_t*>(m_connections[m_worker_id].m_buffer_write);
    buffer_write[0] = message_sz;
    buffer_write[1] = (uint32_t) (force ? RequestType::BATCH_PLAIN_FORCE_YES : RequestType::BATCH_PLAIN_FORCE_NO );
    memcpy(buffer_write + 2, batch, batch_sz * sizeof(library::UpdateInterface::SingleUpdate));
    ssize_t bytes_sent = send(m_connections[m_worker_id].m_fd, (void*) buffer_write, message_sz, /* flags */ 0);
    if(bytes_sent == -1) ERROR_ERRNO("send_request, connection error");
    assert(bytes_sent == message_sz && "Message not fully sent");

    // this is actually more expensive than memcpy into the buffer and doing one send();
//    // send the header
//    assert(m_worker_id >= 0 && m_worker_id < max_num_connections && "Invalid worker id");
//    uint32_t* buffer_write = reinterpret_cast<uint32_t*>(m_connections[m_worker_id].m_buffer_write);
//    buffer_write[0] = header_sz;
//    buffer_write[1] = (uint32_t) RequestType::BATCH_PLAIN;
//    ssize_t bytes_sent = send(m_connections[m_worker_id].m_fd, buffer_write, header_sz, /* flags */ 0);
//    if(bytes_sent == -1) ERROR_ERRNO("send_request, header, connection error");
//    assert(bytes_sent == header_sz && "Message not fully sent");
//
//    // send the body
//    bytes_sent = send(m_connections[m_worker_id].m_fd, (void*) batch, body_sz, /* flags */ 0);
//    if(bytes_sent == -1) ERROR_ERRNO("send_request, body, connection error");
//    assert(bytes_sent == body_sz && "Message not fully sent");

    // wait for the server to process the request
    wait_response();

    if(response()->type() == ResponseType::TIMEOUT){
        assert(force == true && "Only with the flag force == true, the interface is allowed to timeout");
        RAISE_EXCEPTION(library::TimeoutError, "Batch timeout");
    }

    assert(response()->type() == ResponseType::OK);
    return response()->get<bool>(0);
}

void Client::set_timeout(uint64_t seconds){
    const_cast<Client*>(this)->request(RequestType::SET_TIMEOUT, seconds);
    assert(response()->type() == ResponseType::OK);
}

void Client::dump() const {
    const_cast<Client*>(this)->request(RequestType::DUMP_CLIENT);
    assert(response()->type() == ResponseType::OK);
    cout << response()->get_string(0) << endl;
}

void Client::dump(const std::string& path) const {
    const_cast<Client*>(this)->request(RequestType::DUMP_FILE, path);
    assert(response()->type() == ResponseType::OK);
}

void Client::dump_ostream(std::ostream& out) const {
    ERROR("OPERATION NOT SUPPORTED: the output stream is only local");
}

void Client::bfs(uint64_t source_vertex_id, const char* dump2file){
    if(dump2file != nullptr && dump2file[0] == '\0') dump2file = nullptr;


    request(RequestType::BFS, source_vertex_id, dump2file);
    if(response()->type() == ResponseType::NOT_SUPPORTED){
        ERROR("bfs(" << source_vertex_id << ", \"" << dump2file << "\"): operation not supported by the remote interface");
    } else if (response()->type() == ResponseType::ERROR){
        RPC_ERROR(response()->get_string(0));
    } else if (response()->type() == ResponseType::TIMEOUT){
        TIMEOUT_ERROR
    }
    assert(response()->type() == ResponseType::OK);
}

void Client::pagerank(uint64_t num_iterations, double damping_factor, const char* dump2file) {
    request(RequestType::PAGERANK, num_iterations, damping_factor, dump2file);
    if(response()->type() == ResponseType::NOT_SUPPORTED){
        ERROR("pagerank(" << num_iterations << ", " << damping_factor << ", \"" << dump2file << "\"): operation not supported by the remote interface");
    } else if (response()->type() == ResponseType::TIMEOUT){
        TIMEOUT_ERROR
    }
    assert(response()->type() == ResponseType::OK);
}

void Client::wcc(const char* dump2file){
    request(RequestType::WCC, dump2file);
    if(response()->type() == ResponseType::NOT_SUPPORTED){
        ERROR("wcc(\"" << dump2file << "\"): operation not supported by the remote interface");
    } else if (response()->type() == ResponseType::TIMEOUT){
        TIMEOUT_ERROR
    }
    assert(response()->type() == ResponseType::OK);
}

void Client::cdlp(uint64_t max_iterations, const char* dump2file){
    request(RequestType::CDLP, max_iterations, dump2file);

    switch(response()->type()){
    case ResponseType::NOT_SUPPORTED:
        ERROR("cdlp(" << max_iterations << ", \"" << dump2file << "\"): operation not supported by the remote interface");
        break;
    case ResponseType::TIMEOUT:
        TIMEOUT_ERROR
        break;
    case ResponseType::ERROR:
        RPC_ERROR(response()->get_string(0));
        break;
    case ResponseType::OK:
        /* nop */
        break;
    default:
        ERROR("Invalid response type: " << response()->type())
    }
}

void Client::lcc(const char* dump2file){
    request(RequestType::LCC, dump2file);

    switch(response()->type()){
    case ResponseType::NOT_SUPPORTED:
        ERROR("lcc(\"" << dump2file << "\"): operation not supported by the remote interface");
        break;
    case ResponseType::TIMEOUT:
        TIMEOUT_ERROR
        break;
    case ResponseType::ERROR:
        RPC_ERROR(response()->get_string(0));
        break;
    case ResponseType::OK:
        /* nop */
        break;
    default:
        ERROR("Invalid response type: " << response()->type())
    }
}

void Client::sssp(uint64_t source_vertex_id, const char* dump2file){
    request(RequestType::SSSP, source_vertex_id, dump2file);

    switch(response()->type()){
    case ResponseType::NOT_SUPPORTED:
        ERROR("sssp(" << source_vertex_id << ", \"" << dump2file << "\"): operation not supported by the remote interface");
        break;
    case ResponseType::TIMEOUT:
        TIMEOUT_ERROR
        break;
    case ResponseType::ERROR:
        RPC_ERROR(response()->get_string(0));
        break;
    case ResponseType::OK:
        /* nop */
        break;
    default:
        ERROR("Invalid response type: " << response()->type())
    }
}

} // namespace network
