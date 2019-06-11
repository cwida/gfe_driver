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

#include <string>

#include "library/interface.hpp"
#include "message.hpp"

namespace network {

class Client : public library::UpdateInterface, public library::LoaderInterface {
    Client(const Client&) = delete;
    Client& operator=(const Client&) = delete;

    static thread_local int m_worker_id; // keep track which worker
    const std::string m_server_host;
    const int m_server_port;
    static constexpr int max_num_connections = 1024;
    static constexpr size_t buffer_sz = 4096;

    struct ConnectionState {
        int m_fd; // file descriptor for the connection
        char* m_buffer_read; // read buffer
        char* m_buffer_write; // write buffer
    };

    ConnectionState m_connections[max_num_connections]; // keep track of all connections

    /**
     * Open a new connection to the server. Store the file descriptor of the connection in m_fd_connection[m_worker_id]
     */
    void connect();

    /**
     * Close the connection to the server. It affects only the file descriptor referred by m_worker_id
     */
    void disconnect();

    /**
     * Close the connection for the given file descriptor
     */
    void disconnect(int worker_id);

    /**
     * Send the given request to the server
     */
    template<typename... Args>
    void request(RequestType type, Args... args);

    /**
     * Retrieve the current response from the server
     */
    const Response* response() const;

public:
    /**
     * Connect the proxy to the address specified in the configuration
     */
    Client();

    /**
     * Connect the proxy to the server at the given host/port
     */
    Client(const std::string& host, int port);

    /**
     * Destructor
     */
    ~Client();

    /**
     * Cannot dump the content of the graph because the remote prints it in the stdout
     */
    virtual void dump() const override;

    // Proxy to the rest of the functions in the library
    virtual void on_main_init(int num_threads) override;
    virtual void on_thread_init(int thread_id) override;
    virtual void on_thread_destroy(int thread_id) override;
    virtual void on_main_destroy() override;
    virtual uint64_t num_edges() const override;
    virtual uint64_t num_vertices() const override;
    virtual void load(const std::string& path) override;
    virtual bool add_vertex(uint64_t vertex_id) override;
    virtual bool delete_vertex(uint64_t vertex_id) override;
    virtual bool add_edge(graph::WeightedEdge e) override;
    virtual bool delete_edge(graph::Edge e) override;
};

} // namespace network



