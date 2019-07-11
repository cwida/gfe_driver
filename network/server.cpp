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

#include "server.hpp"

#include <arpa/inet.h> // inet_ntoa
#include <cassert>
#include <cstring>
#include <netinet/ip.h> // TCP/IP protocol
#include <sys/select.h>
#include <sys/socket.h>
#include <thread>
#include <unistd.h>

#include "graph/edge.hpp"
#include "library/interface.hpp"
#include "configuration.hpp"
#include "internal.hpp"
#include "message.hpp"

using namespace std;

namespace network {

/*****************************************************************************
 *                                                                           *
 * Server                                                                    *
 *                                                                           *
 *****************************************************************************/

Server::Server(shared_ptr<library::Interface> interface) : Server(interface, configuration().server_port()) { }
Server::Server(shared_ptr<library::Interface> interface, int port) : m_interface(interface), m_port(port){
    LOG("Initialising the server to listen on port: " << port);

    int fd = socket(AF_INET, SOCK_STREAM, 0);
    if(fd < 0) ERROR_ERRNO("Cannot initialise the socket");
    m_server_fd = fd;

    // avoid the error `address already in use'
    bool reuse_value = true;
    setsockopt(m_server_fd, SOL_SOCKET, SO_REUSEADDR, (void*) &reuse_value, sizeof(reuse_value));

    struct sockaddr_in address;
    memset(&address, 0, sizeof(address));
    address.sin_family = AF_INET;
    address.sin_port = m_port;
    address.sin_addr.s_addr = INADDR_ANY; // accept everything
    int rc = bind(m_server_fd, (struct sockaddr*) &address, sizeof(address));
    if(rc != 0) ERROR_ERRNO("Cannot bind the server socket to port: " << m_port);

    // be ready to accept connections
    rc = listen(m_server_fd, SOMAXCONN);
    if(rc != 0) ERROR_ERRNO("Error while attempting to make the socket ready");
}

Server::~Server(){
    // close the file descriptor
    if(m_server_fd >= 0){
        close(m_server_fd);
        m_server_fd = -1;
    }
}

void Server::main_loop(){
    LOG("Server listening to port: " << m_port);

    while(!m_server_stop){
        // Set the timeout to 1 second (it may be changed after each call to select)
        struct timeval timeout;
        timeout.tv_sec = 1;
        timeout.tv_usec = 0;

        // Dummy file descriptor
        fd_set dummy; FD_ZERO(&dummy);

        // Set the file descriptor of the server
        fd_set server_ready;
        FD_ZERO(&server_ready);
        FD_SET(m_server_fd, &server_ready);

        int ready = select(m_server_fd + 1, &server_ready, &dummy, &dummy, &timeout);

        if(ready < 0){
            ERROR_ERRNO("server, select");
        } else if(ready == 0){ // timeout, check the flag m_server_stop
            continue;
        }

        // there is only one fd that can be set, no need to check the bitmask from select
        struct sockaddr_in address;
        socklen_t address_len {0};
        int connection_fd = accept(m_server_fd, (struct sockaddr *) &address, &address_len);
        if(connection_fd < 0)
            ERROR_ERRNO("Cannot establish a connection with a remote client");

        LOG("Connection received from: " << inet_ntoa(address.sin_addr) << ":" << address.sin_port);

        auto t = thread(&ConnectionHandler::execute, new ConnectionHandler(this, connection_fd));
        t.detach(); // do not explicitly wait for the thread to terminate
    }

    LOG("Server terminated");
}


/*****************************************************************************
 *                                                                           *
 * Connection Handler                                                        *
 *                                                                           *
 *****************************************************************************/
Server::ConnectionHandler::ConnectionHandler(Server* instance, int fd) : m_instance(instance), m_fd(fd) { }
Server::ConnectionHandler::~ConnectionHandler() {
    close(m_fd);
    m_fd = -1;
}

void Server::ConnectionHandler::execute(){
    LOG("Remote worker started");

    while(!m_terminate){
        int64_t num_bytes_read { 0 }, recv_bytes { 0 };

        // read the size of the message
        recv_bytes = recv(m_fd, m_buffer_read + num_bytes_read, sizeof(uint32_t), /* flags */ 0);
        if(recv_bytes == -1) ERROR_ERRNO("recv, connection interrupted?");
        assert(recv_bytes == sizeof(uint32_t) && "What the heck have we read?");
        num_bytes_read += recv_bytes;
        int64_t message_sz = static_cast<int64_t>(*(reinterpret_cast<uint32_t*>(m_buffer_read)));
        assert(message_sz + sizeof(uint32_t) < buffer_sz && "Message too long");
        // read the rest of the message
        while(num_bytes_read < message_sz){
            recv_bytes = read(m_fd, m_buffer_read + num_bytes_read, message_sz - num_bytes_read);
            if(recv_bytes == -1) ERROR_ERRNO("recv, connection interrupted?");
            num_bytes_read += recv_bytes;
        }
        assert(num_bytes_read == message_sz && "Message read");

        handle_request();
    }

    LOG("Remote worker terminated");
}

void Server::ConnectionHandler::handle_request(){
    switch(request()->type()){
    case RequestType::TERMINATE_WORKER:
        response(ResponseType::OK);
        m_terminate = true;
        break;
    case RequestType::TERMINATE_SERVER:
        response(ResponseType::OK);
        m_terminate = true;
        m_instance->m_server_stop = true;
        break;
    case RequestType::ON_MAIN_INIT:
        interface()->on_main_init((int) request()->get<int>(0));
        response(ResponseType::OK);
        break;
    case RequestType::ON_THREAD_INIT:
        interface()->on_thread_init((int) request()->get<int>(0));
        response(ResponseType::OK);
        break;
    case RequestType::ON_THREAD_DESTROY:
        interface()->on_thread_destroy((int) request()->get<int>(0));
        response(ResponseType::OK);
        break;
    case RequestType::ON_MAIN_DESTROY:
        interface()->on_main_destroy();
        response(ResponseType::OK);
        break;
    case RequestType::NUM_EDGES: {
        uint64_t num_edges = interface()->num_edges();
        response(ResponseType::OK, num_edges);
    } break;
    case RequestType::NUM_VERTICES: {
        uint64_t num_vertices = interface()->num_vertices();
        response(ResponseType::OK, num_vertices);
    } break;
    case RequestType::HAS_VERTEX: {
        bool result = interface()->has_vertex(request()->get(0));
        response(ResponseType::OK, result);
    } break;
    case RequestType::HAS_EDGE: {
        bool result = interface()->has_edge(request()->get(0), request()->get(1));
        response(ResponseType::OK, result);
    } break;
    case RequestType::GET_WEIGHT: {
        int64_t result = interface()->get_weight(request()->get(0), request()->get(1));
        response(ResponseType::OK, result);
    } break;
    case RequestType::LOAD: {
        library::LoaderInterface* loader_interface = dynamic_cast<library::LoaderInterface*>(interface());
        if(loader_interface == nullptr){
            LOG("Operation not supported by the current interface: " << request()->type());
            response(ResponseType::NOT_SUPPORTED);
        } else {
            string path(request()->buffer());
            LOG("[Server] Attempting to load the graph from path: " << path);
            loader_interface->load(path);
            response(ResponseType::OK);
        }
    } break;
    case RequestType::ADD_VERTEX: {
        library::UpdateInterface* update_interface = dynamic_cast<library::UpdateInterface*>(interface());
        if(update_interface == nullptr){
            LOG("Operation not supported by the current interface: " << request()->type());
            response(ResponseType::NOT_SUPPORTED);
        } else {
            bool result = update_interface->add_vertex(request()->get(0));
            response(ResponseType::OK, result);
        }
    } break;
    case RequestType::DELETE_VERTEX: {
        library::UpdateInterface* update_interface = dynamic_cast<library::UpdateInterface*>(interface());
        if(update_interface == nullptr){
            LOG("Operation not supported by the current interface: " << request()->type());
            response(ResponseType::NOT_SUPPORTED);
        } else {
            bool result = update_interface->delete_vertex(request()->get(0));
            response(ResponseType::OK, result);
        }
    } break;
    case RequestType::ADD_EDGE: {
        library::UpdateInterface* update_interface = dynamic_cast<library::UpdateInterface*>(interface());
        if(update_interface == nullptr){
            LOG("Operation not supported by the current interface: " << request()->type());
            response(ResponseType::NOT_SUPPORTED);
        } else {
            graph::WeightedEdge edge { request()->get(0),  request()->get(1), request()->get<uint32_t>(2)};
            bool result = update_interface->add_edge(edge);
            response(ResponseType::OK, result);
        }
    } break;
    case RequestType::DELETE_EDGE: {
        library::UpdateInterface* update_interface = dynamic_cast<library::UpdateInterface*>(interface());
        if(update_interface == nullptr){
            LOG("Operation not supported by the current interface: " << request()->type());
            response(ResponseType::NOT_SUPPORTED);
        } else {
            graph::Edge edge { request()->get(0),  request()->get(1)};
            bool result = update_interface->delete_edge(edge);
            response(ResponseType::OK, result);
        }
    } break;
    case RequestType::BFS: {
        auto graphalytics = dynamic_cast<library::GraphalyticsInterface*>(interface());
        if(graphalytics == nullptr){
            LOG("Operation not supported by the current interface: " << request()->type());
            response(ResponseType::NOT_SUPPORTED);
        } else {
            string path = request()->get_string(1);
            const char* c_path = path.empty() ? nullptr : path.c_str();
            graphalytics->bfs(request()->get(0), c_path);
            response(ResponseType::OK);
        }
    } break;
    case RequestType::PAGERANK: {
        auto graphalytics = dynamic_cast<library::GraphalyticsInterface*>(interface());
        if(graphalytics == nullptr){
            LOG("Operation not supported by the current interface: " << request()->type());
            response(ResponseType::NOT_SUPPORTED);
        } else {
            string path = request()->get_string(2);
            const char* c_path = path.empty() ? nullptr : path.c_str();
            graphalytics->pagerank(request()->get(0), request()->get<double>(1), c_path);
            response(ResponseType::OK);
        }
    } break;
    case RequestType::WCC: {
        auto graphalytics = dynamic_cast<library::GraphalyticsInterface*>(interface());
        if(graphalytics == nullptr){
            LOG("Operation not supported by the current interface: " << request()->type());
            response(ResponseType::NOT_SUPPORTED);
        } else {
            string path = request()->get_string(0);
            const char* c_path = path.empty() ? nullptr : path.c_str();
            graphalytics->wcc(c_path);
            response(ResponseType::OK);
        }
    } break;
    case RequestType::CDLP: {
        auto graphalytics = dynamic_cast<library::GraphalyticsInterface*>(interface());
        if(graphalytics == nullptr){
            LOG("Operation not supported by the current interface: " << request()->type());
            response(ResponseType::NOT_SUPPORTED);
        } else {
            string path = request()->get_string(1);
            const char* c_path = path.empty() ? nullptr : path.c_str();
            graphalytics->cdlp(request()->get(0), c_path);
            response(ResponseType::OK);
        }
    } break;
    case RequestType::LCC: {
        auto graphalytics = dynamic_cast<library::GraphalyticsInterface*>(interface());
        if(graphalytics == nullptr){
            LOG("Operation not supported by the current interface: " << request()->type());
            response(ResponseType::NOT_SUPPORTED);
        } else {
            string path = request()->get_string(0);
            const char* c_path = path.empty() ? nullptr : path.c_str();
            graphalytics->lcc(c_path);
            response(ResponseType::OK);
        }
    } break;
    case RequestType::SSSP: {
        auto graphalytics = dynamic_cast<library::GraphalyticsInterface*>(interface());
        if(graphalytics == nullptr){
            LOG("Operation not supported by the current interface: " << request()->type());
            response(ResponseType::NOT_SUPPORTED);
        } else {
            string path = request()->get_string(1);
            const char* c_path = path.empty() ? nullptr : path.c_str();
            graphalytics->sssp(request()->get(0), c_path);
            response(ResponseType::OK);
        }
    } break;
    default:
        ERROR("Invalid request type: " << request()->type());
        break;
    }
}


/**
 * Retrieve the request being current processed
 */
const Request* Server::ConnectionHandler::request() const {
    return reinterpret_cast<const Request*>(m_buffer_read);
}

template<typename... Args>
void Server::ConnectionHandler::response(ResponseType type, Args... args){
    new (m_buffer_write) Response(type, std::forward<Args>(args)...);

    uint64_t message_sz = *(reinterpret_cast<uint32_t*>(m_buffer_write));
    assert(message_sz < buffer_sz && "The message is too long");
    ssize_t bytes_sent = send(m_fd, m_buffer_write, message_sz, /* flags */ 0);
    if(bytes_sent == -1) ERROR_ERRNO("send_response, connection error");
    assert(bytes_sent == message_sz && "Message not fully sent");
}

library::Interface* Server::ConnectionHandler::interface(){
    return m_instance->m_interface.get();
}

} // namespace

