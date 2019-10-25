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
#include <unordered_map>
#include <utility>
#include <vector>

#include "common/error.hpp"

class BaseConfiguration; // forward declaration
class ClientConfiguration; // forward declaration
class DriverConfiguration; // forward declaration
class ServerConfiguration; // forward declaration
class StandaloneConfiguration; // forward declaration
namespace common { class Database; } // forward declaration
namespace cxxopts { class Options; } // forward declaration
namespace cxxopts { class ParseResult; } // forward declaration
namespace library { class Interface; } // forward declaration

// Singleton interface
BaseConfiguration& configuration(); // retrieve the current singleton (client, server or standalone)
ClientConfiguration& cfgclient(); // retrieve the singleton for the client configuration
DriverConfiguration& cfgdriver(); // retrieve the singleton for the driver configuration
ServerConfiguration& cfgserver(); // retrieve the singleton for the server configuration
StandaloneConfiguration& cfgstandalone(); // retrieve the singleton for the standalone configuration
void cfgfree(); // invoked at the end to release the configuration

// Generic configuration error
DEFINE_EXCEPTION(ConfigurationError);

// Print the given message to the standard output
extern std::mutex _log_mutex;
#define LOG( msg ) { std::scoped_lock lock(_log_mutex); std::cout << msg << /* flush immediately */ std::endl; }

// Type of counter for the number of threads
enum ThreadsType { THREADS_READ, THREADS_WRITE, THREADS_TOTAL };

/**
 * Base class for the configuration. The actual type can be either ClientConfiguration or ServerConfiguration.
 */
class BaseConfiguration {
    // remove the copy ctors
    BaseConfiguration(const BaseConfiguration& ) = delete;
    BaseConfiguration& operator=(const BaseConfiguration& ) = delete;

    // properties
    common::Database* m_database { nullptr }; // handle to the database
    std::string m_database_path { "" }; // the path where to store the results
    double m_max_weight { 1.0 }; // the maximum weight that can be assigned when reading non weighted graphs
    uint64_t m_seed = 5051789ull; // random seed, used in various places in the experiments

protected:
    // Set the path to the database
    void set_database_path(const std::string& path){ m_database_path = path; }

    // The max weight that can be assigned by graph readers when parsing a non weighted graph
    void set_max_weight(double value);

    // Set the property seed
    void set_seed(uint64_t value){ m_seed = value; }

    // Invoked by save_parameters(); get the list of parameters to store in the database
    // Each subclass should override this method to add its own parameters and invoke the same method in the superclass.
    using param_list_t = std::vector<std::pair<std::string,std::string>>;
    virtual param_list_t list_parameters() const;

    // Add the command line options associated to this class. One definition per subclass.
    virtual void cla_add(cxxopts::Options& options);

    // Parse the command line options associated to the base configuration. One definition per subclass.
    virtual void cla_parse(cxxopts::Options& options, cxxopts::ParseResult& result);

    // The name of the program to show the in the help screen
    virtual std::string cla_name() const;

    // Parse the arguments from the command line
    void parse_command_line_args(int argc, char* argv[]);

public:
    // Constructor
    BaseConfiguration();

    // Destructor
    virtual ~BaseConfiguration();

    // Check whether the configuration/results need to be stored into a database
    bool has_database() const;

    // Retrieve the handle to the database connection, where the final results of the experiments are stored
    common::Database* db();

    // Save the configuration properties into the database
    void save_parameters();

    // Random seed, used in various places in the experiments
    uint64_t seed() const { return m_seed; };

    // Is this a client, the server or standalone program?
    bool is_client() const;
    bool is_server() const;
    bool is_standalone() const;

    // Get the max weight that can be assigned by the reader to
    double max_weight() const { return m_max_weight; }

    // Retrieve the path to the database
    const std::string& get_database_path() const { return m_database_path; }
};

/**
 * The global configuration for the server program (gfe_server).
 * - Initialise (only one time) the singleton instance through ServerConfiguration::initialise(int argc, char* argv[])
 * - Access the singleton instance through the function ::cfgserver();
 * - At the end of the execution, release the singleton with ::cfgfree();
 * The class is not thread safe.
 */
class ServerConfiguration : public BaseConfiguration {
public:
    // remove the copy ctors
    ServerConfiguration(const ServerConfiguration& ) = delete;
    ServerConfiguration& operator=(const ServerConfiguration& ) = delete;
    static constexpr uint32_t DEFAULT_PORT = 18286;

private:
    bool m_graph_directed = true; // whether the graph is undirected or directed
    std::string m_library_name; // the library to test
    std::unique_ptr<library::Interface> (*m_library_factory)(bool directed, uint64_t) {nullptr} ; // function to retrieve an instance of the library `m_library_name'
    uint32_t m_port = DEFAULT_PORT; // the port to wait for new connections

protected:
    // Do not explicitly initialise the configuration, use the method ::initialise();
    ServerConfiguration();

    // Add the command line options associated to this class. One definition per subclass.
    void cla_add(cxxopts::Options& options) override;

    // Parse the command line options associated to the base configuration. One definition per subclass.
    void cla_parse(cxxopts::Options& options, cxxopts::ParseResult& result) override;

    // The name of the program to show the in the help screen
    std::string cla_name() const override;

    // The list of parameters to store into the database
    param_list_t list_parameters() const override;

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

    // Whether the graph is directed or undirected
    bool is_graph_directed() const { return m_graph_directed; }
};

/**
 * Base class for both ClientConfiguration and StandaloneCnfiguration
 */
class DriverConfiguration : public BaseConfiguration {
    // remove the copy ctors
    DriverConfiguration(const DriverConfiguration& ) = delete;
    DriverConfiguration& operator=(const DriverConfiguration& ) = delete;

private:
    uint64_t m_build_frequency { 5 * 60 * 1000 }; // in the aging experiment, the amount of time that must pass before each invocation to #build(), in milliseconds
    double m_coeff_aging { 0.0 }; // coefficient for the additional updates to perform
    double m_ef_vertices = 1; // expansion factor for the vertices in the graph
    double m_ef_edges = 1;  // expansion factor for the edges in the graph
    uint64_t m_num_repetitions { 5 }; // when applicable, how many times the same experiment should be repeated
    int m_num_threads_read { 1 }; // number of threads to use for the read operations
    int m_num_threads_write { 1 }; // number of threads to use for the write (insert/update/delete) operations
    uint64_t m_timeout_seconds { 3600 }; // max time to complete an operation, in seconds (0 => indefinite)
    std::string m_path_graph_to_load; // the file must be accessible to the server

protected:
    void set_build_frequency(uint64_t millisecs);
    void set_coeff_aging(double value); // Set the coefficient for `aging', i.e. how many updates (insertions/deletions) to perform w.r.t. to the size of the loaded graph
    void set_ef_vertices(double value);
    void set_ef_edges(double value);
    void set_num_repetitions(uint64_t value); // Set how many times to repeat the Graphalytics suite of algorithms
    void set_num_thread_read(int value); // Set the number of threads to use in the read operations.
    void set_num_thread_write(int value); // Set the number of threads to use in the write operations.
    void set_timeout(uint64_t seconds); // Set the timeout property
    void set_graph(const std::string& graph); // Set the graph to load and run the experiments

    // Add the command line options associated to this class. One definition per subclass.
    void cla_add(cxxopts::Options& options) override;

    // Parse the command line options associated to the base configuration. One definition per subclass.
    void cla_parse(cxxopts::Options& options, cxxopts::ParseResult& result) override;

    // The list of parameters to store into the database
    param_list_t list_parameters() const override;
public:
    // Constructor
    DriverConfiguration() = default;

    // Coefficient for the surplus of updates to perform (noise) w.r.t. the final graph  to load
    double coefficient_aging() const{ return m_coeff_aging; }

    // Number of repetitions of the same experiment (when applicable)
    uint64_t num_repetitions() const { return m_num_repetitions; }

    // Get the number of threads to use
    int num_threads(ThreadsType type) const;

    // The path for the graph to load
    const std::string& get_path_graph() const { return m_path_graph_to_load; }

    // The budget to complete a Graphalytics algorithm, in seconds (e.g. LCC should terminate by get_timeout_per_operation() seconds)
    uint64_t get_timeout_per_operation() const { return m_timeout_seconds; }

    // Get the expansion factor in the aging experiment for the edges in the graph
    double get_ef_edges() const { return m_ef_edges; }

    // Get the expansion factor in the aging experiment for the vertices in the graph
    double get_ef_vertices() const { return m_ef_vertices; }

    // Get the frequency to build a new snapshot, in milliseconds
    uint64_t get_build_frequency() const{ return m_build_frequency; }
};


/**
 * The global configuration for the client program (gfe_client).
 * - Initialise (only one time) the singleton instance through ClientConfiguration::initialise(int argc, char* argv[])
 * - Access the singleton instance through the function ::cfgclient();
 * - At the end of the execution, release the singleton with ::cfgfree();
 * The class is not thread safe.
 */
class ClientConfiguration : public DriverConfiguration {
    // remove the copy ctors
    ClientConfiguration(const ClientConfiguration& ) = delete;
    ClientConfiguration& operator=(const ClientConfiguration& ) = delete;

private:
    uint64_t m_batch_size = 0; // if > 0, send the updates in batches
    std::string m_experiment = "basic"; // the experiment to execute
    bool m_is_interactive { false }; // interactive while loop
    std::string m_server_host = "localhost";
    uint32_t m_server_port = ServerConfiguration::DEFAULT_PORT;
    bool m_terminate_server_on_exit { false }; // whether to terminate the server after the experiments have been completed

protected:
    // Do not explicitly initialise the configuration, use the method ::initialise();
    ClientConfiguration();

    // Add the command line options associated to this class. One definition per subclass.
    void cla_add(cxxopts::Options& options) override;

    // Parse the command line options associated to the base configuration. One definition per subclass.
    void cla_parse(cxxopts::Options& options, cxxopts::ParseResult& result) override;

    // The name of the program to show the in the help screen
    std::string cla_name() const override;

    // The list of parameters to store into the database
    param_list_t list_parameters() const override;

public:
    // Initialise the server configuration
    static void initialise(int argc, char* argv[]);

    // Destructor
    ~ClientConfiguration();

    // Get the hostname address of the server
    const std::string& get_server_host() const { return m_server_host; }

    // Get the port of the server
    const uint32_t get_server_port() const { return m_server_port; }

    // The remote address for the server, as a pair hostname:port
    std::string get_server_string() const;

    // Interactive mode?
    bool is_interactive() const { return m_is_interactive; }

    // The experiment to execute
    const std::string& get_experiment_name() const { return m_experiment; }

    // Ask whether to shut down the server once the experiment has been performed
    bool is_terminate_server_on_exit() const { return m_terminate_server_on_exit; }

    // Shall we send the updates in batches of the given size?
    uint64_t get_batch_size() const { return m_batch_size; }
};


/**
 * The global configuration for the standalone program (gfe_standalone).
 * - Initialise (only one time) the singleton instance through StandaloneConfiguration::initialise(int argc, char* argv[])
 * - Access the singleton instance through the function ::cfgstandalone();
 * - At the end of the execution, release the singleton with ::cfgfree();
 * The class is not thread safe.
 */
class StandaloneConfiguration : public DriverConfiguration {
    // remove the copy ctors
    StandaloneConfiguration(const ServerConfiguration& ) = delete;
    StandaloneConfiguration& operator=(const StandaloneConfiguration& ) = delete;

    // properties
    bool m_dense_vertices = false; // whether the vertices in the input graph are in the domain [0, N), N = #total vertices
    bool m_graph_directed = true; // whether the graph is undirected or directed
    std::string m_library_name; // the library to test
    std::string m_update_log; // aging experiment through the log file
    std::unique_ptr<library::Interface> (*m_library_factory)(bool directed, uint64_t) {nullptr} ; // function to retrieve an instance of the library `m_library_name'
    bool m_validate_output = false; // whether to validate the execution results of the Graphalytics algorithms

    // Do not explicitly initialise the configuration, use the method ::initialise();
    StandaloneConfiguration();

    // Add the command line options associated to this class. One definition per subclass.
    void cla_add(cxxopts::Options& options) override;

    // Parse the command line options associated to the base configuration. One definition per subclass.
    void cla_parse(cxxopts::Options& options, cxxopts::ParseResult& result) override;

    // The name of the program to show the in the help screen
    std::string cla_name() const override;

    // The list of parameters to store into the database
    param_list_t list_parameters() const override;

public:
    // Initialise the standalone configuration
    static void initialise(int argc, char* argv[]);

    // Destructor
    ~StandaloneConfiguration();

    const std::string& get_library_name() const { return m_library_name; }

    const std::string& get_update_log() const { return m_update_log; }

    // Generate an instance of the graph library to evaluate
    std::unique_ptr<library::Interface> generate_graph_library();

    // Whether the vertices in the input graph are in the domain [0, N), N = #total vertices
    bool has_dense_vertices() const { return m_dense_vertices; }

    // Whether the graph is directed or undirected
    bool is_graph_directed() const { return m_graph_directed; }

    // Whether to validate the execution results of the Graphalytics algorithms
    bool validate_output() const { return m_validate_output; }
};

