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

/*****************************************************************************
 *                                                                           *
 * Initialisation                                                            *
 *                                                                           *
 *****************************************************************************/

thread_local int Client::m_worker_id { 0 };

Client::Client() : Client(configuration().server_host(), configuration().server_port()) {
    ASSERT(configuration().is_remote_client() && "Not expected to provide a remote connection");
}

Client::Client(const std::string& host, int port) : m_server_host(host), m_server_port(port) {
    // reset the content of the connections
    for(int i = 0; i < max_num_connections; i++){
        m_connections[i].m_fd = -1;
        m_connections[i].m_buffer_read = m_connections[i].m_buffer_write = nullptr;
    }

    LOG("Connecting to " << m_server_host << ":" << m_server_port << " ...");
    connect();
    LOG("Connection established!");
}

Client::~Client(){
    // when the client is destroyed, also the server should terminate
    m_worker_id = 0;
    request(RequestType::TERMINATE_SERVER);

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

    m_connections[m_worker_id].m_fd = fd;
    m_connections[m_worker_id].m_buffer_read = (char*) malloc(sizeof(char) * buffer_sz);
    m_connections[m_worker_id].m_buffer_write = (char*) malloc(sizeof(char) * buffer_sz);

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
    uint64_t message_sz = reinterpret_cast<uint32_t*>(buffer)[0];

    // send the request to the server
    ssize_t bytes_sent = send(m_connections[m_worker_id].m_fd, buffer, message_sz, /* flags */ 0);
    if(bytes_sent == -1) ERROR_ERRNO("send_response, connection error");
    assert(bytes_sent == message_sz && "Message not fully sent");

    // receive the reply from the server
    int64_t num_bytes_read { 0 }, recv_bytes { 0 };
    buffer = m_connections[m_worker_id].m_buffer_read;
    recv_bytes = recv(m_connections[m_worker_id].m_fd, buffer + num_bytes_read, sizeof(uint32_t), /* flags */ 0);
    if(recv_bytes == -1) ERROR_ERRNO("recv, connection interrupted?");
    assert(recv_bytes == sizeof(uint32_t) && "What the heck have we read?");
    num_bytes_read += recv_bytes;
    message_sz = *(reinterpret_cast<uint32_t*>(buffer));
    assert(message_sz + sizeof(uint32_t) < buffer_sz && "Message too long");
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
        disconnect();
    }
    m_worker_id = 0;
}

void Client::on_main_destroy() {
    request(RequestType::ON_MAIN_DESTROY);
    assert(response()->type() == ResponseType::OK);
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

bool Client::delete_vertex(uint64_t vertex_id){
    request(RequestType::DELETE_VERTEX, vertex_id);
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

bool Client::delete_edge(graph::Edge e){
    request(RequestType::DELETE_EDGE, e.source(), e.destination());
    if(response()->type() == ResponseType::NOT_SUPPORTED){
        ERROR("delete_edge(" << e.source() << ", " << e.destination() << "): operation not supported by the remote interface");
    }
    assert(response()->type() == ResponseType::OK);
    return response()->get<bool>(0);
}

void Client::dump() const {
    const_cast<Client*>(this)->request(RequestType::DUMP_STDOUT);
    assert(response()->type() == ResponseType::OK);
}

void Client::dump(const std::string& path) const {
    const_cast<Client*>(this)->request(RequestType::DUMP_FILE, path);
    assert(response()->type() == ResponseType::OK);
}

void Client::dump(std::ostream& out) const {
    ERROR("OPERATION NOT SUPPORTED: the output stream is only local");
}

void Client::bfs(uint64_t source_vertex_id, const char* dump2file){
    request(RequestType::BFS, source_vertex_id, dump2file);
    if(response()->type() == ResponseType::NOT_SUPPORTED){
        ERROR("bfs(" << source_vertex_id << ", \"" << dump2file << "\"): operation not supported by the remote interface");
    }
    assert(response()->type() == ResponseType::OK);
}

void Client::pagerank(uint64_t num_iterations, double damping_factor, const char* dump2file) {
    request(RequestType::PAGERANK, num_iterations, damping_factor, dump2file);
    if(response()->type() == ResponseType::NOT_SUPPORTED){
        ERROR("pagerank(" << num_iterations << ", " << damping_factor << ", \"" << dump2file << "\"): operation not supported by the remote interface");
    }
    assert(response()->type() == ResponseType::OK);
}

void Client::wcc(const char* dump2file){
    request(RequestType::WCC, dump2file);
    if(response()->type() == ResponseType::NOT_SUPPORTED){
        ERROR("wcc(\"" << dump2file << "\"): operation not supported by the remote interface");
    }
    assert(response()->type() == ResponseType::OK);
}

void Client::cdlp(uint64_t max_iterations, const char* dump2file){
    request(RequestType::CDLP, max_iterations, dump2file);
    if(response()->type() == ResponseType::NOT_SUPPORTED){
        ERROR("cdlp(" << max_iterations << ", \"" << dump2file << "\"): operation not supported by the remote interface");
    }
    assert(response()->type() == ResponseType::OK);
}

void Client::lcc(const char* dump2file){
    request(RequestType::LCC, dump2file);
    if(response()->type() == ResponseType::NOT_SUPPORTED){
        ERROR("lcc(\"" << dump2file << "\"): operation not supported by the remote interface");
    }
    assert(response()->type() == ResponseType::OK);
}

void Client::sssp(uint64_t source_vertex_id, const char* dump2file){
    request(RequestType::SSSP, source_vertex_id, dump2file);
    if(response()->type() == ResponseType::NOT_SUPPORTED){
        ERROR("sssp(" << source_vertex_id << ", \"" << dump2file << "\"): operation not supported by the remote interface");
    }
    assert(response()->type() == ResponseType::OK);
}

} // namespace network
