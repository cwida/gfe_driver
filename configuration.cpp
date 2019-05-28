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

#include <string>
#include <unistd.h> // sysconf

#include "common/cpu_topology.hpp"
#include "common/database.hpp"
#include "common/quantity.hpp"
#include "common/system.hpp"
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
    m_num_threads = cpu_topology().get_threads(false, false).size();
}
Configuration::~Configuration() {
    delete m_database; m_database = nullptr;
    cout << "Done" << endl;
}

/*****************************************************************************
 *                                                                           *
 *  Command line arguments                                                   *
 *                                                                           *
 *****************************************************************************/
void Configuration::parse_command_line_args(int argc, char* argv[]){
    using namespace cxxopts;

    Options opts(argv[0], "Evaluate the graph libraries");

    opts.add_options("Generic")
        ("d, database", "The path where to store the results", value<string>()->default_value(m_database_path))
        ("e, experiment", "The experiment to execute", value<string>())
        ("h, help", "Show this help menu")
        ("seed", "Random seed used in various places in the experiments", value<uint64_t>()->default_value(to_string(m_seed)))
        ("t, threads", "The number of client threads to use", value<int>()->default_value(to_string(m_num_threads)))
        ("v, verbose", "Print additional messages to the output")
    ;

    try {

        auto result = opts.parse(argc, argv);

        if(result.count("help") > 0){
            cout << opts.help({"Generic"}) << endl;
            exit(0);
        }

        m_database_path = result["database"].as<string>();
        ASSERT( result["threads"].as<int>() >= 0 );
        m_num_threads = result["threads"].as<int>();
        m_seed = result["seed"].as<uint64_t>();
        m_verbose = ( result.count("verbose") > 0 );

    } catch ( argument_incorrect_type& e){
        ERROR(e.what());
    }
}

/************************************************************************\*****
 *                                                                           *
 *  Database                                                                 *
 *                                                                           *
 *****************************************************************************/
common::Database* Configuration::db(){
    if(m_database == nullptr){
        m_database = new Database{m_database_path};
        m_database->create_execution();
    }
    return m_database;
}

void Configuration::save_parameters() {
    using P = pair<string, string>;
    vector<P> params;
    params.push_back(P{"database", m_database_path});
    params.push_back(P{"git_commit", common::git_last_commit()});
    params.push_back(P{"hostname", common::hostname()});
    params.push_back(P{"num_threads", to_string(m_num_threads)});
    params.push_back(P{"seed", to_string(m_seed)});
    params.push_back(P{"verbose", to_string(m_verbose)});

    sort(begin(params), end(params));
    db()->store_parameters(params);
}


