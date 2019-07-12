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

#include "configuration.hpp"

#include <algorithm>
#include <cctype> // tolower
#include <sstream>
#include <string>
#include <unistd.h> // sysconf

#include "common/cpu_topology.hpp"
#include "common/database.hpp"
#include "common/quantity.hpp"
#include "common/system.hpp"
#include "library/interface.hpp"
#include "third-party/cxxopts/cxxopts.hpp"

using namespace common;
using namespace std;

#undef CURRENT_ERROR_TYPE
#define CURRENT_ERROR_TYPE ConfigurationError

/*****************************************************************************
 *                                                                           *
 *  Singleton                                                                *
 *                                                                           *
 *****************************************************************************/
static Configuration singleton;
Configuration& configuration(){ return singleton; }

Configuration::Configuration() {
    m_num_threads_read = m_num_threads_write = cpu_topology().get_threads(false, false).size();
}
Configuration::~Configuration() {
    delete m_database; m_database = nullptr;
}

/*****************************************************************************
 *                                                                           *
 *  Command line arguments                                                   *
 *                                                                           *
 *****************************************************************************/
static string libraries_help_screen(){
    auto libs = library::implementations();
    sort(begin(libs), end(libs), [](const auto& lib1, const auto& lib2){
       return lib1.m_name < lib2.m_name;
    });

    stringstream stream;
    stream << "The library to evaluate: ";
    bool first = true;
    for(const auto& l : libs){
        if (first) first = false; else stream << ", ";
        stream << l.m_name;
    }

    return stream.str();
}

void Configuration::parse_command_line_args(int argc, char* argv[]){
    using namespace cxxopts;

    Options opts(argv[0], "Evaluate the graph libraries");

    opts.add_options("Generic")
        ("a, aging", "The number of additional updates for the aging experiment to perform", value<Quantity>()->default_value(to_string(m_num_aging_updates)))
        ("d, database", "Store the current configuration\result into the given database")
        ("e, experiment", "The experiment to execute", value<string>())
        ("G, graph", "The path to the graph to load", value<string>())
        ("h, help", "Show this help menu")
        ("l, library", libraries_help_screen())
        ("max_weight", "The maximum weight that can be assigned when reading non weighted graphs", value<uint64_t>()->default_value(to_string(m_max_weight)))
        ("r, readers", "The number of client threads to use for the read operations", value<int>()->default_value(to_string(m_num_threads_read)))
        ("seed", "Random seed used in various places in the experiments", value<uint64_t>()->default_value(to_string(m_seed)))
        ("server", "Remote connection, provide a string host:port for the client, and just the port for the server", value<std::string>())
        ("t, threads", "The number of threads to use for both the read and write operations", value<int>()->default_value(to_string(m_num_threads_read + m_num_threads_write)))
        ("v, verbose", "Print additional messages to the output")
        ("w, writers", "The number of client threads to use for the write operations", value<int>()->default_value(to_string(m_num_threads_write)))
    ;

    try {

        auto result = opts.parse(argc, argv);

        if(result.count("help") > 0){
            cout << opts.help({"Generic"}) << endl;
            exit(0);
        }

        m_num_aging_updates = result["aging"].as<Quantity>();
        m_graph_path = result["graph"].as<string>();
        m_seed = result["seed"].as<uint64_t>();
        m_verbose = ( result.count("verbose") > 0 );

        uint64_t max_weight = result["max_weight"].as<uint64_t>();
        if(max_weight <= 0){
            ERROR("Invalid value for the max weight, it must be an integer strictly positive, given: " << max_weight);
        }
        m_max_weight = max_weight;

        // database path
        if( result["database"].count() > 0 ){
            m_database_path = result["database"].as<string>();
        }

        // library to evaluate
        if( result["library"].count() == 0 ){
            ERROR("Missing mandatory argument --library. Which library do you want to evaluate??");
        } else {
            string library_name = result["library"].as<string>();
            transform(begin(library_name), end(library_name), begin(library_name), ::tolower); // make it lower case
            auto libs = library::implementations();
            auto library_found = find_if(begin(libs), end(libs), [&library_name](const auto& candidate){
                return library_name == candidate.m_name;
            });
            if(library_found == end(libs)){ ERROR("Library not recognised: `" << result["library"].as<string>() << "'"); }
            m_library_name = library_found->m_name;
            m_library_factory = library_found->m_factory;
        }

        // number of threads
        if( result["threads"].count() > 0) {
           ASSERT( result["threads"].as<int>() >= 0 );
           m_num_threads_read = m_num_threads_write = result["threads"].as<int>();
        }
        if( result["readers"].count() > 0) {
            ASSERT( result["readers"].as<int>() >= 0 );
            m_num_threads_read = result["readers"].as<int>();
        }
        if( result["writers"].count() > 0 ){
            ASSERT( result["writers"].as<int>() >= 0 );
            m_num_threads_write = result["writers"].as<int>();
        }

        // network connection
        if( result["server"].count() > 0 ){
            string param_server = result["server"].as<string>();
            auto pos_colon = param_server.find(':');
            string param_port = param_server;

            if(pos_colon != string::npos){
                m_server_host = param_server.substr(0, pos_colon);
                param_port = param_server.substr(pos_colon +1);
            }

            // parse the port number
            try {
                int port = stoi(param_port); // throws invalid_argument if the conversion cannot be performed
                if(port <= 0 || port >= (1>>16)){ throw invalid_argument("invalid port number"); }
                m_server_port = port;
            } catch (invalid_argument& e){
                ERROR("Invalid parameter --server: " << param_server << ". The port number cannot be recognised: `" << param_port << "'");
            }
        };

    } catch ( argument_incorrect_type& e){
        ERROR(e.what());
    }
}

/*****************************************************************************
 *                                                                           *
 *  Database                                                                 *
 *                                                                           *
 *****************************************************************************/
bool Configuration::has_database() const {
    return !m_database_path.empty();
}

common::Database* Configuration::db(){
    if(m_database == nullptr && has_database()){
        m_database = new Database{m_database_path};
        m_database->create_execution();
    }
    return m_database;
}

void Configuration::save_parameters() {
    if(db() == nullptr) ERROR("Path where to store the results not set");

    using P = pair<string, string>;
    vector<P> params;
    params.push_back(P{"database", m_database_path});
    params.push_back(P{"git_commit", common::git_last_commit()});
    params.push_back(P{"graph", graph()});
    params.push_back(P{"hostname", common::hostname()});
    params.push_back(P{"library", m_library_name});
    params.push_back(P{"num_aging_updates", to_string(m_num_aging_updates)});
    params.push_back(P{"num_threads_read", to_string(m_num_threads_read)});
    params.push_back(P{"num_threads_write", to_string(m_num_threads_write)});
    params.push_back(P{"seed", to_string(m_seed)});
    params.push_back(P{"verbose", to_string(m_verbose)});

    // remote connection
    if(is_remote_client()){
        params.push_back(P{"remote_role", "client"});
        params.push_back(P{"server_host", server_host()});
        params.push_back(P{"server_port", to_string(server_port())});
    } else if (is_remote_server()){
        params.push_back(P{"remote_role", "server"});
        params.push_back(P{"server_port", to_string(server_port())});
    }

    sort(begin(params), end(params));
    db()->store_parameters(params);
}

/*****************************************************************************
 *                                                                           *
 *  Properties                                                               *
 *                                                                           *
 *****************************************************************************/
int Configuration::num_threads(ThreadsType type) const {
    switch(type){
    case THREADS_READ:
        return m_num_threads_read;
    case THREADS_WRITE:
        return m_num_threads_write;
    case THREADS_TOTAL:
        return m_num_threads_read + m_num_threads_write;
    default:
        ERROR("Invalid thread type: " << ((int) type));
    }
}

bool Configuration::is_remote_client() const {
    return !m_server_host.empty() && m_server_port > 0;
}

bool Configuration::is_remote_server() const {
    return m_server_host.empty() && m_server_port > 0;
}

std::string Configuration::server_host() const {
    return m_server_host;
}

int Configuration::server_port() const {
    return m_server_port;
}

std::unique_ptr<library::Interface> Configuration::generate_graph_library() {
    return m_library_factory();
}

