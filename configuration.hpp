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
#include <memory>
#include <mutex>

#include "common/error.hpp"

class Configuration; // forward declaration
class ClientConfiguration; // forward declaration
class ServerConfiguration; // forward declaration
namespace common { class Database; } // forward declaration
namespace library { class Interface; } // forward declaration

// Singleton interface
Configuration& configuration();
ClientConfiguration& cfgclient();
ServerConfiguration& cfgserver();
void cfgfree(); // invoked at the end to release the configuration

// Generic configuration error
DEFINE_EXCEPTION(ConfigurationError);

// Print the given message only if we are in verbose mode
extern std::mutex _log_mutex;

#define LOG( msg ) { std::scoped_lock lock(_log_mutex); std::cout << msg << /* flush immediately */ std::endl; }

// Type of counter for the number of threads
enum ThreadsType { THREADS_READ, THREADS_WRITE, THREADS_TOTAL };

/**
 * Base class for the configuration. The actual type can be either ClientConfiguration or ServerConfiguration.
 */
class Configuration {
    // remove the copy ctors
    Configuration(const Configuration& ) = delete;
    Configuration& operator=(const Configuration& ) = delete;

    // properties
    common::Database* m_database { nullptr }; // handle to the database
    std::string m_database_path { "" }; // the path where to store the results
    double m_max_weight { 1024.0 }; // the maximum weight that can be assigned when reading non weighted graphs
    uint64_t m_seed = 5051789ull; // random seed, used in various places in the experiments

protected:
    // Set the path to the database
    void set_database_path(const std::string& path){ m_database_path = path; }

    // The max weight that can be assigned by graph readers when parsing a non weighted graph
    void set_max_weight(double value);

    // Set the property seed
    void set_seed(uint64_t value){ m_seed = value; }

public:
    // Constructor
    Configuration();

    // Destructor
    virtual ~Configuration();

    // Check whether the configuration/results need to be stored into a database
    bool has_database() const;

    // Retrieve the handle to the database connection, where the final results of the experiments are stored
    common::Database* db();

    // Save the configuration properties into the database
    virtual void save_parameters(); // nop for this class, expected to be overriden by the subclasses

    // Random seed, used in various places in the experiments
    uint64_t seed() const { return m_seed; };

    // Is this a client or the server?
    bool is_client() const;
    bool is_server() const;

    // Get the max weight that can be assigned by the reader to
    double max_weight() const { return m_max_weight; }

    // Retrieve the path to the database
    const std::string& get_database_path() const { return m_database_path; }
};

class ServerConfiguration : public Configuration {
public:
    // remove the copy ctors
    ServerConfiguration(const ServerConfiguration& ) = delete;
    ServerConfiguration& operator=(const ServerConfiguration& ) = delete;
    static constexpr uint32_t DEFAULT_PORT = 18286;

private:
    bool m_graph_directed = true; // whether the graph is undirected or directed
    std::string m_library_name; // the library to test
    std::unique_ptr<library::Interface> (*m_library_factory)(bool directed) {nullptr} ; // function to retrieve an instance of the library `m_library_name'
    uint32_t m_port = DEFAULT_PORT; // the port to wait for new connections

protected:
    // Do not explicitly initialise the configuration, use the method ::initialise();
    ServerConfiguration();

    // Parse the arguments from the command line
    void parse_command_line_args(int argc, char* argv[]);

public:
    // Initialise the server configuration
    static void initialise(int argc, char* argv[]);

    // Destructor
    ~ServerConfiguration();

    // Get the port to wait for new client connections
    uint32_t get_port() const { return m_port; }

    const std::string& get_library_name() const { return m_library_name; }

    // Generate an instance of the graph library to evaluate
    std::unique_ptr<library::Interface> generate_graph_library();

    // Save the configuration properties into the database
    virtual void save_parameters() override;

    // Whether the graph is directed or undirected
    bool is_graph_directed() const { return m_graph_directed; }
};


class ClientConfiguration : public Configuration {
    // remove the copy ctors
    ClientConfiguration(const ClientConfiguration& ) = delete;
    ClientConfiguration& operator=(const ClientConfiguration& ) = delete;

private:
    uint64_t m_batch_size = 0; // if > 0, send the updates in batches
    double m_coeff_aging { 0.0 }; // coefficient for the additional updates to perform
    double m_ef_vertices = 1; // expansion factor for the vertices in the graph
    double m_ef_edges = 1;  // expansion factor for the edges in the graph
    std::string m_experiment; // the experiment to execute
    bool m_is_interactive { false }; // interactive while loop
    uint64_t m_num_repetitions { 5 }; // when applicable, how many times the same experiment should be repeated
    int m_num_threads_read { 1 }; // number of threads to use for the read operations
    int m_num_threads_write { 1 }; // number of threads to use for the write (insert/update/delete) operations
    uint64_t m_timeout_seconds { 3600 }; // max time to complete an operation, in seconds (0 => indefinite)
    std::string m_path_graph_to_load; // the file must be accessible to the server
    std::string m_server_host = "localhost";
    uint32_t m_server_port = ServerConfiguration::DEFAULT_PORT;
    bool m_terminate_server_on_exit { false }; // whether to terminate the server after the experiments have been completed

protected:
    // Do not explicitly initialise the configuration, use the method ::initialise();
    ClientConfiguration();

    // Parse the arguments from the command line
    void parse_command_line_args(int argc, char* argv[]);

public:
    // Initialise the server configuration
    static void initialise(int argc, char* argv[]);

    // Destructor
    ~ClientConfiguration();

    // Save the configuration properties into the database
    virtual void save_parameters() override;

    // Get the hostname address of the server
    const std::string& get_server_host() const { return m_server_host; }

    // Get the port of the server
    const uint32_t get_server_port() const { return m_server_port; }

    // The remote address for the server, as a pair hostname:port
    std::string get_server_string() const;

    // Interactive mode?
    bool is_interactive() const { return m_is_interactive; }

    // Get the number of threads to use
    int num_threads(ThreadsType type) const;

    // Coefficient for the surplus of updates to perform (noise) w.r.t. the final graph  to load
    double coefficient_aging() const{ return m_coeff_aging; }

    // Number of repetitions of the same experiment (when applicable)
    uint64_t num_repetitions() const { return m_num_repetitions; }

    // The path for the graph to load
    const std::string& get_path_graph() const { return m_path_graph_to_load; }

    // The experiment to execute
    const std::string& get_experiment_name() const { return m_experiment; }

    bool is_terminate_server_on_exit() const { return m_terminate_server_on_exit; }

    uint64_t get_timeout_per_operation() const { return m_timeout_seconds; }

    // Shall we send the updates in batches of the given size?
    uint64_t get_batch_size() const { return m_batch_size; }

    // Get the expansion factor in the aging experiment for the edges in the graph
    double get_ef_edges() const { return m_ef_edges; }

    // Get the expansion factor in the aging experiment for the vertices in the graph
    double get_ef_vertices() const { return m_ef_vertices; }
};


