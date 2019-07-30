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
 *  Base class                                                               *
 *                                                                           *
 *****************************************************************************/
mutex _log_mutex;

static Configuration* singleton {nullptr};
Configuration& configuration(){
    if(singleton == nullptr){ ERROR("Configuration not initialised"); }
    return *singleton;
}

void cfgfree(){ delete singleton; singleton = nullptr; }

Configuration::Configuration() {
//    m_num_threads_read = m_num_threads_write = cpu_topology().get_threads(false, false).size();
}

Configuration::~Configuration() {
    delete m_database; m_database = nullptr;
}

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

void Configuration::set_max_weight(double value){
    if(value <= 0) ERROR("Invalid value for max weight: " << value << ". Expected a positive value");
    m_max_weight = value;
}

bool Configuration::is_client() const {
    return dynamic_cast<const ClientConfiguration*>(this) != nullptr;
}

bool Configuration::is_server() const {
    return dynamic_cast<const ServerConfiguration*>(this) != nullptr;
}


/*****************************************************************************
 *                                                                           *
 *  Client configuration                                                     *
 *                                                                           *
 *****************************************************************************/
void ClientConfiguration::initialise(int argc, char* argv[]){
    if(singleton != nullptr){ ERROR("Global configuration already initialised!"); }
    auto cfg = new ClientConfiguration();
    singleton = cfg;
    cfg->parse_command_line_args(argc, argv);
}
ClientConfiguration::ClientConfiguration() {  }
ClientConfiguration::~ClientConfiguration() {  }

ClientConfiguration& cfgclient(){
    auto cfg = dynamic_cast<ClientConfiguration*>(singleton);
    if(cfg == nullptr){ ERROR("Configuration not initialised or not in client mode"); }
    return *cfg;
}

void ClientConfiguration::parse_command_line_args(int argc, char* argv[]){
    using namespace cxxopts;

    Options opts(argv[0], "Client program for the GFE");

    opts.add_options("Generic")
        ("a, aging", "The number of additional updates for the aging experiment to perform", value<Quantity>()->default_value(to_string(m_num_aging_updates)))
        ("c, connect", "The server address, in the form hostname:port. The default is " + get_server_string(), value<string>())
        ("d, database", "Store the current configuration value into the a sqlite3 database at the given location", value<string>())
        ("e, experiment", "The experiment to execute", value<string>())
        ("h, help", "Show this help menu")
        ("G, graph", "The path to the graph to load", value<string>())
        ("i, interactive", "Show the command loop to interact with the server")
        ("max_weight", "The maximum weight that can be assigned when reading non weighted graphs", value<double>()->default_value(to_string(max_weight())))
        ("R, repetitions", "The number of repetitions of the same experiment (where applicable)", value<uint64_t>()->default_value(to_string(m_num_repetitions)))
        ("r, readers", "The number of client threads to use for the read operations", value<int>()->default_value(to_string(m_num_threads_read)))
        ("seed", "Random seed used in various places in the experiments", value<uint64_t>()->default_value(to_string(seed())))
        ("t, threads", "The number of threads to use for both the read and write operations", value<int>()->default_value(to_string(m_num_threads_read + m_num_threads_write)))
//        ("v, verbose", "Print additional messages to the output")
        ("w, writers", "The number of client threads to use for the write operations", value<int>()->default_value(to_string(m_num_threads_write)))
    ;

    try {
        auto result = opts.parse(argc, argv);

        if(result.count("help") > 0){
            cout << opts.help({"Generic"}) << endl;
            exit(EXIT_SUCCESS);
        }

        set_seed( result["seed"].as<uint64_t>() );
//        set_verbose( ( result.count("verbose") > 0 ) );
        set_max_weight( result["max_weight"].as<double>() );
        if( result["database"].count() > 0 ){ set_database_path( result["database"].as<string>() ); }

        if( result["connect"].count() > 0 ){
            string param_server = result["connect"].as<string>();
            auto pos_colon = param_server.find(':');
            string param_port = param_server;

            if(pos_colon != string::npos){
                m_server_host = param_server.substr(0, pos_colon);
                param_port = param_server.substr(pos_colon +1);
            }

            // parse the port number
            try {
                int port = stoi(param_port); // throws invalid_argument if the conversion cannot be performed
                if(port <= 0 || port >= (1<<16)){ throw invalid_argument("invalid port number"); }
                m_server_port = port;
            } catch (invalid_argument& e){
                ERROR("Invalid parameter --connect: " << param_server << ". The port number cannot be recognised: `" << param_port << "'");
            }
        };

        m_is_interactive = (result["interactive"].count() > 0);
        m_num_aging_updates = result["aging"].as<Quantity>();

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

        // the graph to work with
        if( result["graph"].count() > 0 ){
            m_path_graph_to_load = result["graph"].as<string>();
        }

        // the experiment to execute
        if( result["experiment"].count() > 0 ){
            m_experiment = result["experiment"].as<string>();
            transform(begin(m_experiment), end(m_experiment), begin(m_experiment), ::tolower);
        }

        uint64_t num_repetitions = result["repetitions"].as<uint64_t>();
        if(num_repetitions <= 0) ERROR("Invalid value for the parameter --repetitions: " << num_repetitions << ". Expected a positive value");
        m_num_repetitions = num_repetitions;
    } catch ( argument_incorrect_type& e){
        ERROR(e.what());
    }
}

void ClientConfiguration::save_parameters() {
    if(db() == nullptr) ERROR("Path where to store the results not set");

    using P = pair<string, string>;
    vector<P> params;
    params.push_back(P{"aging", to_string(m_num_aging_updates)});
    params.push_back(P{"database", get_database_path()});
    if(!get_experiment_name().empty()) params.push_back(P{"experiment", get_experiment_name()});
    params.push_back(P{"git_commit", common::git_last_commit()});
    if(!get_path_graph().empty()){ params.push_back(P{"graph", get_path_graph()}); }
    params.push_back(P{"interactive", to_string( is_interactive() )});
    params.push_back(P{"hostname", common::hostname()});
    params.push_back(P{"num_aging_updates", to_string(m_num_aging_updates)});
    params.push_back(P{"num_repetitions", to_string(m_num_repetitions)});
    params.push_back(P{"num_threads_read", to_string(m_num_threads_read)});
    params.push_back(P{"num_threads_write", to_string(m_num_threads_write)});
    params.push_back(P{"role", "client"});
    params.push_back(P{"seed", to_string(seed())});
    params.push_back(P{"server_host", get_server_host()});
    params.push_back(P{"server_port", to_string(get_server_port())});
//    params.push_back(P{"verbose", to_string(verbose())});


    sort(begin(params), end(params));
    db()->store_parameters(params);
}

string ClientConfiguration::get_server_string() const {
    stringstream ss;
    ss << get_server_host() << ":" << get_server_port();
    return ss.str();
}

int ClientConfiguration::num_threads(ThreadsType type) const {
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

/*****************************************************************************
 *                                                                           *
 *  Server configuration                                                     *
 *                                                                           *
 *****************************************************************************/
void ServerConfiguration::initialise(int argc, char* argv[]){
    if(singleton != nullptr){ ERROR("Global configuration already initialised!"); }
    auto cfg = new ServerConfiguration();
    singleton = cfg;
    cfg->parse_command_line_args(argc, argv);
}
ServerConfiguration::ServerConfiguration() {  }
ServerConfiguration::~ServerConfiguration() {  }

ServerConfiguration& cfgserver(){
    auto cfg = dynamic_cast<ServerConfiguration*>(singleton);
    if(cfg == nullptr){ ERROR("Configuration not initialised or not in server mode"); }
    return *cfg;
}

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

void ServerConfiguration::parse_command_line_args(int argc, char* argv[]){
    using namespace cxxopts;

    Options opts(argv[0], "Server program for the GFE");

    opts.add_options("Generic")
        ("d, database", "Store the current configuration value into the a sqlite3 database at the given location", value<string>())
        ("h, help", "Show this help menu")
        ("l, library", libraries_help_screen(), value<string>())
        ("max_weight", "The maximum weight that can be assigned when reading non weighted graphs", value<double>()->default_value(to_string(max_weight())))
        ("p, port", "The port where to accept remote connections from the clients", value<uint32_t>()->default_value( to_string(get_port()) ))
        ("seed", "Random seed used in various places in the experiments", value<uint64_t>()->default_value(to_string(seed())))
        ("u, undirected", "Is the graph undirected? By default, it's considered directed.")
//        ("v, verbose", "Print additional messages to the output")
    ;

    try {

        auto result = opts.parse(argc, argv);

        if(result.count("help") > 0){
            cout << opts.help({"Generic"}) << endl;
            exit(EXIT_SUCCESS);
        }

        set_seed( result["seed"].as<uint64_t>() );
//        set_verbose( ( result.count("verbose") > 0 ) );
        set_max_weight( result["max_weight"].as<double>() );
        if( result["database"].count() > 0 ){ set_database_path( result["database"].as<string>() ); }

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

        m_port = result["port"].as<uint32_t>();

        if( result["undirected"].count() > 0 ){
            m_graph_directed = false;
        }

    } catch ( argument_incorrect_type& e){
        ERROR(e.what());
    }
}

void ServerConfiguration::save_parameters() {
    if(db() == nullptr) ERROR("Path where to store the results not set");

    using P = pair<string, string>;
    vector<P> params;
    params.push_back(P{"database", get_database_path()});
    params.push_back(P{"directed", to_string(is_graph_directed())});
    params.push_back(P{"git_commit", common::git_last_commit()});
    params.push_back(P{"hostname", common::hostname()});
    params.push_back(P{"library", m_library_name});
    params.push_back(P{"port", to_string(get_port())});
    params.push_back(P{"role", "server"});
    params.push_back(P{"seed", to_string(seed())});
//    params.push_back(P{"verbose", to_string(verbose())});

    sort(begin(params), end(params));
    db()->store_parameters(params);
}

std::unique_ptr<library::Interface> ServerConfiguration::generate_graph_library() {
    return m_library_factory(is_graph_directed());
}

///*****************************************************************************
// *                                                                           *
// *  Command line arguments                                                   *
// *                                                                           *
// *****************************************************************************/
//static string libraries_help_screen(){
//    auto libs = library::implementations();
//    sort(begin(libs), end(libs), [](const auto& lib1, const auto& lib2){
//       return lib1.m_name < lib2.m_name;
//    });
//
//    stringstream stream;
//    stream << "The library to evaluate: ";
//    bool first = true;
//    for(const auto& l : libs){
//        if (first) first = false; else stream << ", ";
//        stream << l.m_name;
//    }
//
//    return stream.str();
//}
//
//void Configuration::parse_command_line_args(int argc, char* argv[]){
//    using namespace cxxopts;
//
//    Options opts(argv[0], "Evaluate the graph libraries");
//
//    opts.add_options("Generic")
//        ("a, aging", "The number of additional updates for the aging experiment to perform", value<Quantity>()->default_value(to_string(m_num_aging_updates)))
//        ("d, database", "Store the current configuration\result into the given database")
//        ("e, experiment", "The experiment to execute", value<string>())
//        ("G, graph", "The path to the graph to load", value<string>())
//        ("h, help", "Show this help menu")
//        ("l, library", libraries_help_screen())
//        ("max_weight", "The maximum weight that can be assigned when reading non weighted graphs", value<uint64_t>()->default_value(to_string(m_max_weight)))
//        ("r, readers", "The number of client threads to use for the read operations", value<int>()->default_value(to_string(m_num_threads_read)))
//        ("seed", "Random seed used in various places in the experiments", value<uint64_t>()->default_value(to_string(m_seed)))
//        ("server", "Remote connection, provide a string host:port for the client, and just the port for the server", value<std::string>())
//        ("t, threads", "The number of threads to use for both the read and write operations", value<int>()->default_value(to_string(m_num_threads_read + m_num_threads_write)))
//        ("v, verbose", "Print additional messages to the output")
//        ("w, writers", "The number of client threads to use for the write operations", value<int>()->default_value(to_string(m_num_threads_write)))
//    ;
//    opts.a
//
//    try {
//
//        auto result = opts.parse(argc, argv);
//
//        if(result.count("help") > 0){
//            cout << opts.help({"Generic"}) << endl;
//            exit(0);
//        }
//
//
//        m_graph_path = result["graph"].as<string>();
//        m_seed = result["seed"].as<uint64_t>();
//        m_verbose = ( result.count("verbose") > 0 );
//
//        uint64_t max_weight = result["max_weight"].as<uint64_t>();
//        if(max_weight <= 0){
//            ERROR("Invalid value for the max weight, it must be an integer strictly positive, given: " << max_weight);
//        }
//        m_max_weight = max_weight;
//
//        // database path
//        if( result["database"].count() > 0 ){
//            m_database_path = result["database"].as<string>();
//        }
//
//        // library to evaluate
//        if( result["library"].count() == 0 ){
//            ERROR("Missing mandatory argument --library. Which library do you want to evaluate??");
//        } else {
//            string library_name = result["library"].as<string>();
//            transform(begin(library_name), end(library_name), begin(library_name), ::tolower); // make it lower case
//            auto libs = library::implementations();
//            auto library_found = find_if(begin(libs), end(libs), [&library_name](const auto& candidate){
//                return library_name == candidate.m_name;
//            });
//            if(library_found == end(libs)){ ERROR("Library not recognised: `" << result["library"].as<string>() << "'"); }
//            m_library_name = library_found->m_name;
//            m_library_factory = library_found->m_factory;
//        }
//
//        // number of threads
//        if( result["threads"].count() > 0) {
//           ASSERT( result["threads"].as<int>() >= 0 );
//           m_num_threads_read = m_num_threads_write = result["threads"].as<int>();
//        }
//        if( result["readers"].count() > 0) {
//            ASSERT( result["readers"].as<int>() >= 0 );
//            m_num_threads_read = result["readers"].as<int>();
//        }
//        if( result["writers"].count() > 0 ){
//            ASSERT( result["writers"].as<int>() >= 0 );
//            m_num_threads_write = result["writers"].as<int>();
//        }
//
//        // network connection
//        if( result["server"].count() > 0 ){
//            string param_server = result["server"].as<string>();
//            auto pos_colon = param_server.find(':');
//            string param_port = param_server;
//
//            if(pos_colon != string::npos){
//                m_server_host = param_server.substr(0, pos_colon);
//                param_port = param_server.substr(pos_colon +1);
//            }
//
//            // parse the port number
//            try {
//                int port = stoi(param_port); // throws invalid_argument if the conversion cannot be performed
//                if(port <= 0 || port >= (1>>16)){ throw invalid_argument("invalid port number"); }
//                m_server_port = port;
//            } catch (invalid_argument& e){
//                ERROR("Invalid parameter --server: " << param_server << ". The port number cannot be recognised: `" << param_port << "'");
//            }
//        };
//
//    } catch ( argument_incorrect_type& e){
//        ERROR(e.what());
//    }
//}

