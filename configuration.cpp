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
#include <random>
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
 *  Helpers                                                                  *
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

/*****************************************************************************
 *                                                                           *
 *  Base class                                                               *
 *                                                                           *
 *****************************************************************************/
mutex _log_mutex;

static BaseConfiguration* singleton {nullptr};
BaseConfiguration& configuration(){
    if(singleton == nullptr){
//        LOG("[configuration] initialising a test configuration ... ");
        singleton = new BaseConfiguration();
    }
    return *singleton;
}

void cfgfree(){ delete singleton; singleton = nullptr; }

BaseConfiguration::BaseConfiguration() {
//    m_num_threads_read = m_num_threads_write = cpu_topology().get_threads(false, false).size();
}

BaseConfiguration::~BaseConfiguration() {
    delete m_database; m_database = nullptr;
}

void BaseConfiguration::parse_command_line_args(int argc, char* argv[]){
    using namespace cxxopts;

    Options opts(argv[0], cla_name());
    cla_add(opts);

    try {
        auto result = opts.parse(argc, argv);
        cla_parse(opts, result);

    } catch ( argument_incorrect_type& e){
        ERROR(e.what());
    }
}

void BaseConfiguration::cla_add(cxxopts::Options& options){
    using namespace cxxopts;

    options.add_options("Generic")
        ("d, database", "Store the current configuration value into the a sqlite3 database at the given location", value<string>())
        ("h, help", "Show this help menu")
        ("max_weight", "The maximum weight that can be assigned when reading non weighted graphs", value<double>()->default_value(to_string(max_weight())))
        ("seed", "Random seed used in various places in the experiments", value<uint64_t>()->default_value(to_string(seed())))
    ;
}

void BaseConfiguration::cla_parse(cxxopts::Options& options, cxxopts::ParseResult& result){
    if(result.count("help") > 0){
        cout << options.help({"Generic"}) << endl;
        exit(EXIT_SUCCESS);
    }

    set_seed( result["seed"].as<uint64_t>() );
    set_max_weight( result["max_weight"].as<double>() );
    if( result["database"].count() > 0 ){ set_database_path( result["database"].as<string>() ); }
}

string BaseConfiguration::cla_name() const { // this ought to be specialised in the subclasseses, e.g. "GFE Server", "GFE Client", etc.
    return "GFE";
}

bool BaseConfiguration::has_database() const {
    return !m_database_path.empty();
}

common::Database* BaseConfiguration::db(){
    if(m_database == nullptr && has_database()){
        m_database = new Database{m_database_path};
        auto params = m_database->create_execution();
        // random value with no semantic, the aim is to simplify the work of ./automerge.pl
        // in recognising duplicate entries
        params.add("magic", (uint64_t) std::random_device{}());
    }
    return m_database;
}

void BaseConfiguration::set_max_weight(double value){
    if(value <= 0) ERROR("Invalid value for max weight: " << value << ". Expected a positive value");
    m_max_weight = value;
}

bool BaseConfiguration::is_client() const {
    return dynamic_cast<const ClientConfiguration*>(this) != nullptr;
}

bool BaseConfiguration::is_server() const {
    return dynamic_cast<const ServerConfiguration*>(this) != nullptr;
}

bool BaseConfiguration::is_standalone() const {
    return dynamic_cast<const StandaloneConfiguration*>(this) != nullptr;
}

void BaseConfiguration::save_parameters() {
    if(db() == nullptr) ERROR("Path where to store the results not set");
    auto params = list_parameters();
    sort(begin(params), end(params));
    db()->store_parameters(params);
}

auto BaseConfiguration::list_parameters() const -> param_list_t {
    using P = pair<string, string>;

    param_list_t params;
    params.push_back(P{"database", get_database_path()});
    params.push_back(P{"git_commit", common::git_last_commit()});
    params.push_back(P{"hostname", common::hostname()});
    params.push_back(P{"seed", to_string(seed())});
    return params;
}


/*****************************************************************************
 *                                                                           *
 *  Client/standalone common base class                                      *
 *                                                                           *
 *****************************************************************************/
DriverConfiguration& cfgdriver(){
    auto cfg = dynamic_cast<DriverConfiguration*>(singleton);
    if(cfg == nullptr){ ERROR("Configuration not initialised or not in driver mode"); }
    return *cfg;
}

void DriverConfiguration::cla_add(cxxopts::Options& options){
    using namespace cxxopts;

    BaseConfiguration::cla_add(options);

    options.add_options("Generic")
        ("a, aging", "The number of additional updates for the aging experiment to perform", value<double>()->default_value("0"))
        ("efe", "Expansion factor for the edges in the graph", value<double>()->default_value(to_string(get_ef_edges())))
        ("efv", "Expansion factor for the vertices in the graph", value<double>()->default_value(to_string(get_ef_vertices())))
        ("G, graph", "The path to the graph to load", value<string>())
        ("R, repetitions", "The number of repetitions of the same experiment (where applicable)", value<uint64_t>()->default_value(to_string(num_repetitions())))
        ("r, readers", "The number of client threads to use for the read operations", value<int>()->default_value(to_string(num_threads(THREADS_READ))))
        ("timeout", "Set the maximum time for an operation to complete, in seconds", value<uint64_t>()->default_value(to_string(get_timeout_per_operation())))
        ("t, threads", "The number of threads to use for both the read and write operations", value<int>()->default_value(to_string(num_threads(THREADS_TOTAL))))
        ("w, writers", "The number of client threads to use for the write operations", value<int>()->default_value(to_string(num_threads(THREADS_WRITE))))
    ;
}

void DriverConfiguration::cla_parse(cxxopts::Options& options, cxxopts::ParseResult& result){
    BaseConfiguration::cla_parse(options, result);

    // number of threads
    if( result["threads"].count() > 0) {
        int value = result["threads"].as<int>();

        set_num_thread_write(value);
    }
    if( result["readers"].count() > 0) {
        set_num_thread_read( result["readers"].as<int>() );
    }
    if( result["writers"].count() > 0 ){
        set_num_thread_write( result["writers"].as<int>() );
    }

    // the graph to work with
    if( result["graph"].count() > 0 ){
        set_graph( result["graph"].as<string>() );
    }

    set_coeff_aging( result["aging"].as<double>() );
    set_num_repetitions( result["repetitions"].as<uint64_t>() );
    set_timeout( result["timeout"].as<uint64_t>() );
    set_ef_vertices( result["efe"].as<double>() );
    set_ef_vertices( result["efv"].as<double>() );
}

auto DriverConfiguration::list_parameters() const -> param_list_t {
    using P = pair<string, string>;

    param_list_t params = BaseConfiguration::list_parameters();
    params.push_back(P{"aging", to_string(coefficient_aging())});
    params.push_back(P{"ef_edges", to_string(get_ef_edges())});
    params.push_back(P{"ef_vertices", to_string(get_ef_vertices())});
    if(!get_path_graph().empty()){ params.push_back(P{"graph", get_path_graph()}); }
    params.push_back(P{"num_repetitions", to_string(num_repetitions())});
    params.push_back(P{"num_threads_read", to_string(num_threads(ThreadsType::THREADS_READ))});
    params.push_back(P{"num_threads_write", to_string(num_threads(ThreadsType::THREADS_WRITE))});
    params.push_back(P{"timeout", to_string(get_timeout_per_operation())});
    return params;
}

void DriverConfiguration::set_coeff_aging(double value){
    if(value < 0 || (value > 0 && value < 1)){
        ERROR("The parameter aging is invalid. It must be >= 1.0: " << value);
    }
    m_coeff_aging = value;
}

void DriverConfiguration::set_ef_vertices(double value){
    m_ef_vertices = value;
}

void DriverConfiguration::set_ef_edges(double value){
    m_ef_edges = value;
}

void DriverConfiguration::set_num_repetitions(uint64_t value) {
    m_num_repetitions = value; // accept 0 as value
}

void DriverConfiguration::set_num_thread_read(int value){
    ASSERT( value >= 0 );
    m_num_threads_read = value;
}

void DriverConfiguration::set_num_thread_write(int value){
    ASSERT( value >= 0 );
    m_num_threads_write = value;
}

void DriverConfiguration::set_timeout(uint64_t seconds) {
    m_timeout_seconds = seconds;
}

void DriverConfiguration::set_graph(const std::string& graph){
    m_path_graph_to_load = graph;
}

int DriverConfiguration::num_threads(ThreadsType type) const {
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

void ClientConfiguration::cla_add(cxxopts::Options& opts){
    using namespace cxxopts;

    DriverConfiguration::cla_add(opts);

    opts.add_options("Generic")
        ("b, batch", "Send the updates in batches of the given size", value<Quantity>())
        ("c, connect", "The server address, in the form hostname:port. The default is " + get_server_string(), value<string>())
        ("e, experiment", "The experiment to execute", value<string>()->default_value(m_experiment))
        ("i, interactive", "Show the command loop to interact with the server")
        ("p, port", "Specify the port of the remote server", value<uint32_t>())
        ("terminate_server_on_exit", "Terminate the server after the client finished")
    ;
}

void ClientConfiguration::cla_parse(cxxopts::Options& options, cxxopts::ParseResult& result) {
    using namespace cxxopts;

    DriverConfiguration::cla_parse(options, result);

    bool port_specified_in_connect = false;
    if( result["connect"].count() > 0 ){
        string param_server = result["connect"].as<string>();
        auto pos_colon = param_server.find(':');
        if(pos_colon == string::npos) { // port not specified
            m_server_host = param_server;
        } else {
            m_server_host = param_server.substr(0, pos_colon);
            string param_port = param_server.substr(pos_colon +1);

            // parse the port number
            try {
                int port = stoi(param_port); // throws invalid_argument if the conversion cannot be performed
                if(port <= 0 || port >= (1<<16)){ throw invalid_argument("invalid port number"); }
                m_server_port = port;
                port_specified_in_connect = true;
            } catch (invalid_argument& e){
                ERROR("Invalid parameter --connect: " << param_server << ". The port number cannot be recognised: `" << param_port << "'");
            }
        }
    };
    if( result["port"].count() > 0 ){ // alias for -c <port>
        // parse the port number
        uint32_t param_port = result["port"].as<uint32_t>();
        if(port_specified_in_connect && param_port != m_server_port){
            ERROR("The parameter port has been specified twice with different values -c " << m_server_port << " and -p " << param_port);
        } else if (param_port <= 0 || param_port >= (1<<16)){
            throw invalid_argument("invalid port number");
        } else {
            m_server_port = param_port;
        }
    }

    m_is_interactive = (result["interactive"].count() > 0);

    // the experiment to execute
    if( result["experiment"].count() > 0 ){
        m_experiment = result["experiment"].as<string>();
        transform(begin(m_experiment), end(m_experiment), begin(m_experiment), ::tolower);
    }

    m_terminate_server_on_exit = result["terminate_server_on_exit"].count() > 0;

    if(result["batch"].count() > 0){
        m_batch_size = result["batch"].as<Quantity>();
    }
}

string ClientConfiguration::cla_name() const {
    return "Client program for the GFE";
}

auto ClientConfiguration::list_parameters() const -> param_list_t {
    using P = pair<string, string>;
    auto params = DriverConfiguration::list_parameters();
    params.push_back(P{"batch", to_string(m_batch_size)});
    if(!get_experiment_name().empty()) params.push_back(P{"experiment", get_experiment_name()});
    params.push_back(P{"interactive", to_string( is_interactive() )});
    params.push_back(P{"role", "client"});
    params.push_back(P{"server_host", get_server_host()});
    params.push_back(P{"server_port", to_string(get_server_port())});
    params.push_back(P{"terminate_server_on_exit", to_string(is_terminate_server_on_exit())});
    return params;
}

string ClientConfiguration::get_server_string() const {
    stringstream ss;
    ss << get_server_host() << ":" << get_server_port();
    return ss.str();
}

/*****************************************************************************
 *                                                                           *
 *  Standalone configuration                                                 *
 *                                                                           *
 *****************************************************************************/
void StandaloneConfiguration::initialise(int argc, char* argv[]){
    if(singleton != nullptr){ ERROR("Global configuration already initialised!"); }
    auto cfg = new StandaloneConfiguration();
    singleton = cfg;
    cfg->parse_command_line_args(argc, argv);
}
StandaloneConfiguration::StandaloneConfiguration() {  }
StandaloneConfiguration::~StandaloneConfiguration() {  }

StandaloneConfiguration& cfgstandalone(){
    auto cfg = dynamic_cast<StandaloneConfiguration*>(singleton);
    if(cfg == nullptr){ ERROR("Configuration not initialised or not in standalone mode"); }
    return *cfg;
}

void StandaloneConfiguration::cla_add(cxxopts::Options& options) {
    using namespace cxxopts;

    DriverConfiguration::cla_add(options);

    options.add_options("Generic")
        ("l, library", libraries_help_screen(), value<string>())
        ("u, undirected", "Is the graph undirected? By default, it's considered directed.")
    ;
}

void StandaloneConfiguration::cla_parse(cxxopts::Options& options, cxxopts::ParseResult& result) {
    DriverConfiguration::cla_parse(options, result);

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

    if( result["undirected"].count() > 0 ){
        m_graph_directed = false;
    }
}

string StandaloneConfiguration::cla_name() const {
    return "GFE driver";
}

auto StandaloneConfiguration::list_parameters() const -> param_list_t {
    using P = pair<string, string>;

    auto params = DriverConfiguration::list_parameters();
    params.push_back(P{"directed", to_string(is_graph_directed())});
    params.push_back(P{"library", get_library_name()});
    params.push_back(P{"role", "standalone"});

    return params;
}

std::unique_ptr<library::Interface> StandaloneConfiguration::generate_graph_library() {
    return m_library_factory(is_graph_directed());
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

void ServerConfiguration::cla_add(cxxopts::Options& options) {
    using namespace cxxopts;

    BaseConfiguration::cla_add(options);

    options.add_options("Generic")
        ("l, library", libraries_help_screen(), value<string>())
        ("p, port", "The port where to accept remote connections from the clients", value<uint32_t>()->default_value( to_string(get_port()) ))
        ("u, undirected", "Is the graph undirected? By default, it's considered directed.")
    ;
}

void ServerConfiguration::cla_parse(cxxopts::Options& options, cxxopts::ParseResult& result) {
    BaseConfiguration::cla_parse(options, result);

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
}

string ServerConfiguration::cla_name() const {
    return "Server program for the GFE";
}

auto ServerConfiguration::list_parameters() const -> param_list_t {
    using P = pair<string, string>;

    auto params = BaseConfiguration::list_parameters();
    params.push_back(P{"directed", to_string(is_graph_directed())});
    params.push_back(P{"library", m_library_name});
    params.push_back(P{"port", to_string(get_port())});
    params.push_back(P{"role", "server"});

    return params;
}

std::unique_ptr<library::Interface> ServerConfiguration::generate_graph_library() {
    return m_library_factory(is_graph_directed());
}

