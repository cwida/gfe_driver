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

#include <cstdlib>
#include <iostream>

#include "common/error.hpp"
#include "library/interface.hpp"
#include "network/server.hpp"

#include "configuration.hpp"

using namespace std;

static void run_server(int argc, char* argv[]){
    cout << "[server] Init configuration ... " << endl;
    ServerConfiguration::initialise(argc, argv);

    if(configuration().has_database()){
        cout << "[server] Save the current configuration properties in " << configuration().get_database_path() << endl;
        configuration().save_parameters();
    }

    cout << "[server] Init the graph library \"" << cfgserver().get_library_name() << "\", directed graph: " << boolalpha << cfgserver().is_graph_directed() << endl;
    std::shared_ptr<library::Interface> impl2eval { cfgserver().generate_graph_library() };

    cout << "[server] Init che connection handler ..." << endl;
    network::Server server { impl2eval };
    server.handle_signals();

    cout << "[server] Start the event loop ..." << endl;
    server.main_loop();

    cout << "[server] Done" << endl;
}

int main(int argc, char* argv[]){
    int rc = 0;
    try {
        run_server(argc, argv);
    } catch(common::Error& e){
        cerr << e << endl;
        cerr << "Server terminating due to exception..." << endl;
        rc = 1;
    }

    cfgfree(); // release the configuration object

    return rc;
}


