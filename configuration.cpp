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
#include <cmath>
#include <cstdlib>
#include <random>
#include <sstream>
#include <string>
#include <unistd.h> // sysconf

#include "common/cpu_topology.hpp"
#include "common/database.hpp"
#include "common/filesystem.hpp"
#include "common/quantity.hpp"
#include "common/system.hpp"
#include "experiment/graphalytics.hpp"
#include "library/interface.hpp"
#include "reader/graphlog_reader.hpp"
#include "third-party/cxxopts/cxxopts.hpp"

using namespace common;
using namespace std;

#undef CURRENT_ERROR_TYPE
#define CURRENT_ERROR_TYPE ::gfe::ConfigurationError

namespace gfe {

/*****************************************************************************
 *                                                                           *
 *  Helpers                                                                  *
 *                                                                           *
 *****************************************************************************/
mutex _log_mutex;

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
 *  Singleton                                                                *
 *                                                                           *
 *****************************************************************************/
static Configuration g_configuration;
Configuration& configuration(){ return g_configuration; }

/*****************************************************************************
 *                                                                           *
 *  Initialisation                                                           *
 *                                                                           *
 *****************************************************************************/
Configuration::Configuration(){

}

Configuration::~Configuration(){
    delete m_database; m_database = nullptr;
}

void Configuration::initialiase(int argc, char* argv[]){
    using namespace cxxopts;

    Options options(argv[0], "GFE Driver");

    options.add_options("Generic")
        ("a, aging", "The number of additional updates for the aging experiment to perform", value<double>()->default_value("0"))
        ("aging_step_size", "The step of each recording for the measured progress in the Aging2 experiment. Valid values are 0.1, 0.25, 0.5 and 1.0", value<double>()->default_value("1"))
        ("aging_timeout", "Force terminating the aging experiment after four hours")
        ("blacklist", "Comma separated list of graph algorithms to blacklist and do not execute", value<string>())
        ("build_frequency", "The frequency to build a new snapshot in the aging experiment (default: disabled)", value<DurationQuantity>())
        ("d, database", "Store the current configuration value into the a sqlite3 database at the given location", value<string>())
        ("efe", "Expansion factor for the edges in the graph", value<double>()->default_value(to_string(get_ef_edges())))
        ("efv", "Expansion factor for the vertices in the graph", value<double>()->default_value(to_string(get_ef_vertices())))
        ("G, graph", "The path to the graph to load", value<string>())
        ("h, help", "Show this help menu")
        ("latency", "Measure the latency of inserts/updates, report the average, median, std. dev. and 90/95/97/99 percentiles")
        ("l, library", libraries_help_screen(), value<string>())
        ("log", "Repeat the log of updates specified in the given file", value<string>())
        ("max_weight", "The maximum weight that can be assigned when reading non weighted graphs", value<double>()->default_value(to_string(max_weight())))
        ("omp", "Maximum number of threads that can be used by OpenMP (0 = do not change)", value<int>()->default_value(to_string(num_threads_omp())))
        ("R, repetitions", "The number of repetitions of the same experiment (where applicable)", value<uint64_t>()->default_value(to_string(num_repetitions())))
        ("r, readers", "The number of client threads to use for the read operations", value<int>()->default_value(to_string(num_threads(THREADS_READ))))
        ("seed", "Random seed used in various places in the experiments", value<uint64_t>()->default_value(to_string(seed())))
        ("t, threads", "The number of threads to use for both the read and write operations", value<int>()->default_value(to_string(num_threads(THREADS_TOTAL))))
        ("timeout", "Set the maximum time for an operation to complete, in seconds", value<uint64_t>()->default_value(to_string(get_timeout_per_operation())))
        ("u, undirected", "Is the graph undirected? By default, it's considered directed.")
        ("v, validate", "Whether to validate the output results of the Graphalytics algorithms")
        ("w, writers", "The number of client threads to use for the write operations", value<int>()->default_value(to_string(num_threads(THREADS_WRITE))))
    ;

    try {
        auto result = options.parse(argc, argv);

        if(result.count("help") > 0){
            cout << options.help({"Generic"}) << endl;
            ::exit(EXIT_SUCCESS);
        }

        if(result["log"].count() > 0){
            m_update_log = result["log"].as<string>();
            if(!common::filesystem::exists(m_update_log)){ ERROR("Option --log \"" << m_update_log << "\", the file does not exist"); }

            // verify that the properties from the log file are not also specified via command line
            if(result["aging"].count() > 0) { ERROR("Cannot specify the option --aging together with the log file"); }
            if(result["efe"].count() > 0) { ERROR("Cannot specify the option --efe together with the log file"); }
            if(result["efv"].count() > 0) { ERROR("Cannot specify the option --efv together with the log file"); }
            if(result["max_weight"].count() > 0) { ERROR("Cannot specify the option --max_weight together with the log file"); }

            // read the properties from the log file
            auto log_properties = reader::graphlog::parse_properties(m_update_log);
            if(log_properties.find("aging_coeff") == log_properties.end()) { ERROR("The log file `" << m_update_log << "' does not contain the expected property 'aging_coeff'"); }
            set_coeff_aging(stod(log_properties["aging_coeff"]));
            if(log_properties.find("ef_edges") == log_properties.end()) { ERROR("The log file `" << m_update_log << "' does not contain the expected property 'ef_edges'"); }
            set_ef_edges(stod(log_properties["ef_edges"]));
            if(log_properties.find("ef_vertices") == log_properties.end()) { ERROR("The log file `" << m_update_log << "' does not contain the expected property 'ef_vertices'"); }
            set_ef_vertices(stod(log_properties["ef_vertices"]));
            if(log_properties.find("max_weight") != log_properties.end()){ // otherwise assume the default
                set_max_weight(stod(log_properties["max_weight"]));
            }

            // validate that the graph in the input log file matches the same graph specified in the command line
            auto it_path_graph = log_properties.find("input_graph");
            if(it_path_graph != log_properties.end() && result["graph"].count() > 0){
                auto name_log = common::filesystem::filename(it_path_graph->second);
                auto name_param = common::filesystem::filename(common::filesystem::absolute_path(result["graph"].as<string>()));
                if(name_log != name_param){
                    ERROR("The log file is based on the graph `" << name_log << "', while the parameter --graph refers to `" << name_param << "'");
                }
            }
        }

        set_seed( result["seed"].as<uint64_t>() );

        if( result["max_weight"].count() > 0)
            set_max_weight( result["max_weight"].as<double>() );

        if( result["database"].count() > 0 ){ set_database_path( result["database"].as<string>() ); }

        // number of threads
        if( result["threads"].count() > 0) {
            int value = result["threads"].as<int>();

            set_num_threads_write(value);
        }
        if( result["readers"].count() > 0) {
            set_num_threads_read( result["readers"].as<int>() );
        }
        if( result["writers"].count() > 0 ){
            set_num_threads_write( result["writers"].as<int>() );
        }

        // the graph to work with
        if( result["graph"].count() > 0 ){
            set_graph( result["graph"].as<string>() );
        }

        if(result["aging"].count() > 0)
            set_coeff_aging( result["aging"].as<double>() );

        set_num_repetitions( result["repetitions"].as<uint64_t>() );
        set_timeout( result["timeout"].as<uint64_t>() );

        if( result["efe"].count() > 0 )
            set_ef_edges( result["efe"].as<double>() );

        if( result["efv"].count() > 0 )
            set_ef_vertices( result["efv"].as<double>() );

        if( result["build_frequency"].count() > 0 ){
            set_build_frequency( result["build_frequency"].as<DurationQuantity>().as<chrono::milliseconds>().count() );
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

        if( result["undirected"].count() > 0 ){
            m_graph_directed = false;
        }

        if( result["validate"].count() > 0 ){
            m_validate_inserts = true;
            m_validate_output = true;
        }

        if ( result["omp"].count() > 0 ){
            set_num_threads_omp( result["omp"].as<int>() );
        }

        if ( result["aging_step_size"].count() > 0 ){
            set_aging_step_size( result["aging_step_size"].as<double>() );
        }

        if( result["blacklist"].count() > 0 ){
            string algorithm;
            stringstream ss(result["blacklist"].as<string>());
            while(getline(ss, algorithm, ',')){
                // trim
                algorithm.erase(algorithm.begin(), std::find_if(algorithm.begin(), algorithm.end(), [](int ch) { return !std::isspace(ch); }));
                algorithm.erase(std::find_if(algorithm.rbegin(), algorithm.rend(), [](int ch) { return !std::isspace(ch); }).base(), algorithm.end());

                // to lower case
                std::transform(algorithm.begin(), algorithm.end(), algorithm.begin(),[](unsigned char c){ return std::tolower(c); });

                m_blacklist.push_back(algorithm);
            }

            std::sort(m_blacklist.begin(), m_blacklist.end());
        } // blacklist

        m_measure_latency = result["latency"].count() > 0;
        m_timeout_aging2 = result["aging_timeout"].count() > 0;
    } catch ( argument_incorrect_type& e){
        ERROR(e.what());
    }
}

/*****************************************************************************
 *                                                                           *
 *  Properties                                                               *
 *                                                                           *
 *****************************************************************************/

bool Configuration::has_database() const {
    return !m_database_path.empty();
}

::common::Database* Configuration::db(){
    if(m_database == nullptr && has_database()){
        m_database = new Database{m_database_path};
        auto params = m_database->create_execution();
        // random value with no semantic, the aim is to simplify the work of ./automerge.pl
        // in recognising duplicate entries
        params.add("magic", (uint64_t) std::random_device{}());
    }
    return m_database;
}

void Configuration::set_max_weight(double value){
    if(value <= 0) ERROR("Invalid value for max weight: " << value << ". Expected a positive value");
    m_max_weight = value;
}


void Configuration::set_coeff_aging(double value){
    if(value < 0 || (value > 0 && value < 1)){
        ERROR("The parameter aging is invalid. It must be >= 1.0: " << value);
    }
    m_coeff_aging = value;
}

void Configuration::set_ef_vertices(double value){
    m_ef_vertices = value;
}

void Configuration::set_ef_edges(double value){
    m_ef_edges = value;
}

void Configuration::set_num_repetitions(uint64_t value) {
    m_num_repetitions = value; // accept 0 as value
}

void Configuration::set_num_threads_omp(int value){
    ASSERT( value >= 0 );
#if !defined(HAVE_OPENMP)
    if(value > 0) ERROR("Cannot set the maximum number of threads to use with OpenMP: the driver was not configured with support of OpenMP");
#endif

    m_num_threads_omp = value;
}

void Configuration::set_num_threads_read(int value){
    ASSERT( value >= 0 );
#if !defined(HAVE_OPENMP)
    if(value > 0) ERROR("Cannot set the maximum number of threads to use: the driver was not configured with support of OpenMP");
#endif

    m_num_threads_read = value;
}

void Configuration::set_num_threads_write(int value){
    ASSERT( value > 0 );
    m_num_threads_write = value;
}

void Configuration::set_timeout(uint64_t seconds) {
    m_timeout_seconds = seconds;
}

void Configuration::set_graph(const std::string& graph){
    m_path_graph_to_load = graph;
}

void Configuration::set_build_frequency( uint64_t millisecs ){
    m_build_frequency = millisecs;
}

void Configuration::set_aging_step_size( double value ){
    if(value <= 0 || value > 1){
        ERROR("Invalid value for the aging step size. It must be in (0, 1]. Value given: " << value);
    }
    if ( ::ceil(1.0/value) != ::floor(1.0/value) ){
        ERROR("Value for the step size currently not supported: " << value << ". Expected a value such as 1/(step size) is an integer");
    }
    m_step_size_recordings = value;
}

uint64_t Configuration::get_num_recordings_per_ops() const {
    double step_size = get_aging_step_size();
    if ( ::ceil(1.0/step_size) != ::floor(1.0/step_size) ){
        ERROR("Value for the step size currently not supported: " << step_size << ". Expected a value such as 1/(step size) is an integer");
    }
    return 1.0/step_size;
}

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

int Configuration::num_threads_omp() const {
    return m_num_threads_omp;
}

bool Configuration::is_aging2_timeout_set() const {
    return m_timeout_aging2;
}

std::unique_ptr<library::Interface> Configuration::generate_graph_library() {
    return m_library_factory(is_graph_directed());
}

void Configuration::do_blacklist(bool& property_enabled, const char* name) const {
    if(property_enabled){
        if(find(begin(m_blacklist), end(m_blacklist), name) != end(m_blacklist)){
            string name_upper_case = name;
            std::transform(name_upper_case.begin(), name_upper_case.end(), name_upper_case.begin(),[](unsigned char c){ return std::toupper(c); });
            LOG("> Ignore " << name_upper_case << ", algorithm blacklisted");
            property_enabled = false;
        }
    }
}

void Configuration::blacklist(gfe::experiment::GraphalyticsAlgorithms& algorithms) const {
    do_blacklist(algorithms.bfs.m_enabled, "bfs");
    do_blacklist(algorithms.cdlp.m_enabled, "cdlp");
    do_blacklist(algorithms.lcc.m_enabled, "lcc");
    do_blacklist(algorithms.pagerank.m_enabled, "pagerank");
    do_blacklist(algorithms.sssp.m_enabled, "sssp");
    do_blacklist(algorithms.wcc.m_enabled, "wcc");
}

/*****************************************************************************
 *                                                                           *
 *  Save parameters                                                          *
 *                                                                           *
 *****************************************************************************/

void Configuration::save_parameters() {
    if(db() == nullptr) ERROR("Path where to store the results not set");

    using P = pair<string, string>;
    std::vector<P> params;
    params.push_back(P{"database", get_database_path()});
    params.push_back(P{"git_commit", common::git_last_commit()});
    params.push_back(P{"hostname", common::hostname()});
    params.push_back(P("max_weight", to_string(max_weight())));
    params.push_back(P{"seed", to_string(seed())});
    params.push_back(P{"aging", to_string(coefficient_aging())});
    params.push_back(P{"aging_step_size", to_string(get_aging_step_size())});
    params.push_back(P{"aging_timeout", to_string(is_aging2_timeout_set())});
    params.push_back(P{"build_frequency", to_string(get_build_frequency())}); // milliseconds
    params.push_back(P{"ef_edges", to_string(get_ef_edges())});
    params.push_back(P{"ef_vertices", to_string(get_ef_vertices())});
    if(!get_path_graph().empty()){ params.push_back(P{"graph", get_path_graph()}); }
    params.push_back(P{"measure_latency", to_string(measure_latency())});
    params.push_back(P{"num_repetitions", to_string(num_repetitions())});
    params.push_back(P{"num_threads_omp", to_string(num_threads_omp())});
    params.push_back(P{"num_threads_read", to_string(num_threads(ThreadsType::THREADS_READ))});
    params.push_back(P{"num_threads_write", to_string(num_threads(ThreadsType::THREADS_WRITE))});
    params.push_back(P{"timeout", to_string(get_timeout_per_operation())});
    params.push_back(P{"directed", to_string(is_graph_directed())});
    params.push_back(P{"library", get_library_name()});
    if(!get_update_log().empty()) {
        // version 1: uniform distribution
        // version 2: log file, follow the same node degree distribution of the input graph
        // version 3: use #add_edge_v2
        params.push_back(P{"aging_impl", "version_3"});
        params.push_back(P{"log", get_update_log()});
    }
    params.push_back(P{"role", "standalone"});
    params.push_back(P{"validate_inserts", to_string(validate_inserts())});
    params.push_back(P{"validate_output", to_string(validate_output())});

    if(!m_blacklist.empty()){
        stringstream ss;
        for(auto& s: m_blacklist){
            if(ss.tellp() > 0){ ss << ", "; }
            ss << s;
        }
        params.push_back(P{"blacklist", ss.str()});
    }

    sort(begin(params), end(params));
    db()->store_parameters(params);
}

} // namespace
