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
#include <iostream>

#include "common/error.hpp"

class Configuration; // forward declaration
namespace common { class Database; } // forward declaration

// Singleton interface
Configuration& configuration();

// Generic configuration error
DEFINE_EXCEPTION(ConfigurationError);

// Print the given message only if we are in verbose mode
#define LOG( msg ) { if(configuration().verbose()){ std::cout << msg << "\n"; } }

// Type of counter for the number of threads
enum ThreadsType { THREADS_READ, THREADS_WRITE, THREADS_TOTAL };

class Configuration {
    // remove the copy ctors
    Configuration(const Configuration& ) = delete;
    Configuration& operator=(const Configuration& ) = delete;


    // properties
    common::Database* m_database { nullptr }; // handle to the database
    std::string m_database_path { "results.sqlite3" }; // the path where to store the results
    std::string m_graph_path { "" };
    uint64_t m_num_aging_updates { 0 }; // number of additional updates to perform
    int m_num_threads_read { 0 }; // number of threads to use for the read operations
    int m_num_threads_write { 0 }; // number of threads to use for the write (insert/update/delete) operations
    uint64_t m_seed = 5051789ull; // random seed, used in various places in the experiments
    std::string m_server_host = ""; // the hostname for the remote server
    int m_server_port = -1; // the port of the remote server
    bool m_verbose { false }; // verbose mode?

public:
    // Do not explicitly initialise the configuration, use the method ::configuration();
    Configuration();

    // Destructor
    ~Configuration();

    // Parse the command line arguments
    void parse_command_line_args(int argc, char* argv[]);

    // Retrieve the handle to the database connection, where the final results of the experiments are stored
    common::Database* db();

    // The path of the graph to load
    std::string graph() const { return m_graph_path; }

    // Save the configuration into the database
    void save_parameters();

    // Random seed, used in various places in the experiments
    uint64_t seed() const { return m_seed; };

    // Check whether we are in verbose mode, to print additional message to the output
    bool verbose() const { return m_verbose; }

    // Get the number of threads to use
    int num_threads(ThreadsType type) const;

    // Do we conduct the experiment on a remote server?
    bool is_remote_client() const;

    // Check whether this process is a remote server
    bool is_remote_server() const;

    // Get the hostname of the remote host
    std::string server_host() const;

    // Retrieve the port of the remote server
    int server_port() const;
};



