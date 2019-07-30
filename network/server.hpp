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

#include <atomic>
#include <memory>

#include "message.hpp"

namespace library { class Interface; } // forward decl.

namespace network {

class Server {
    std::shared_ptr<library::Interface> m_interface; // the interface we are serving
    const int m_port; // server port
    int m_server_fd {-1}; // file descriptor used by the server to listen for connections
    bool m_server_stop { false }; // flag to stop the server accepting connections
    std::atomic<int> m_num_active_connections = 0;

    class ConnectionHandler {
        Server* m_instance;
        int m_fd;
        static constexpr size_t buffer_sz = 4096; // capacity of the internal buffers
        char m_buffer_read[buffer_sz]; // read buffer (for requests)
        char m_buffer_write[buffer_sz]; // write buffer (for responses)
        bool m_terminate { false }; // flag to signal to terminate the handler

        /**
         * Send the given response to the client
         */
        template<typename... Args>
        void response(ResponseType type, Args... args);

        /**
         * Send the given message to the client
         */
        void send_message(char* raw_message);

        /**
         * Retrieve the request being current processed
         */
        const Request* request() const;

        /**
         * Process a single request
         */
        void handle_request();

        /**
         * The library we are evaluating
         */
        library::Interface* interface();

    public:
        ConnectionHandler(Server* instance, int fd);

        /**
         * Destructor
         */
        ~ConnectionHandler();

        /**
         * Handle all remote requests from the associated file descriptor
         */
        void execute();
    };
    friend class ConnectionHandler;

public:
    Server(std::shared_ptr<library::Interface> interface);
    Server(std::shared_ptr<library::Interface> interface, int port);

    /**
     * Destructor
     */
    ~Server();

    /**
     * Execute the server
     */
    void main_loop();
};

}


