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

int main(int argc, char* argv[]){
    // Parse the command line arguments
    try {
        configuration().parse_command_line_args(argc, argv);
    } catch (common::Error& e){
        cerr << "Error raised while parsing the command line arguments:\n";
        cerr << e << "\n";
        exit(EXIT_FAILURE);
    }

    // Save the experiment's configuration into the output SQLite3 database
    try {
        configuration().save_parameters();
    } catch (common::Error& e){
        cerr << "Error raised while saving the configuration into the database:\n";
        cerr << e << "\n";
        exit(EXIT_FAILURE);
    }

    // Init the library to evaluate
    std::shared_ptr<library::Interface> interface;
    try {
        interface = configuration().generate_graph_library();
    } catch (common::Error& e){
        cerr << "Error raise while generate the library interface:\n";
        cerr << e << "\n";
        exit(EXIT_FAILURE);
    }

    // Remote server mode ?
    if (configuration().is_remote_server()){
        try {
            network::Server server{interface};
            server.main_loop();

        } catch(common::Error& e){
            cerr << "Server error:\n";
            cerr << e << "\n";
            exit(EXIT_FAILURE);
        }

        return 0; // exit
    }

    return 0;
}

