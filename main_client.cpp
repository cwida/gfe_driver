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
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <iostream>

#include "common/error.hpp"
#include "library/interface.hpp"
#include "network/client.hpp"

#include "configuration.hpp"

using namespace std;

/*****************************************************************************
 *                                                                           *
 *  DEBUG                                                                    *
 *                                                                           *
 *****************************************************************************/
#define DEBUG
#define COUT_DEBUG_FORCE(msg) { std::cout << "[main_client::" << __FUNCTION__ << "] " << msg << std::endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif


/*****************************************************************************
 *                                                                           *
 *  Interactive mode                                                         *
 *                                                                           *
 *****************************************************************************/
namespace {
// forward declaration
class Argument;
class Command;
static Command parse_command();

DEFINE_EXCEPTION(ParseCommandError);
#define PARSE_COMMAND_ERROR(message) RAISE_EXCEPTION(ParseCommandError, message)

class Argument {
    friend Command parse_command();
    const string m_argument;

    Argument(const char* start, int length) : m_argument(start, length) {}

public:
    string to_string() const{ return m_argument; }
    operator string() const { return to_string(); }
};

class Command {
    friend Command parse_command();

    string m_command;
    vector<Argument> m_arguments;

    Command(){};
public:
    const std::string& get_command() const{ return m_command; }

    std::string to_string() const {
        stringstream ss;
        ss << m_command << "(";
        bool first = true;
        for(auto& a: m_arguments){
            if(first) { first = false; } else { ss << ", "; }
            ss << a.to_string();
        }
        ss << ")";
        return ss.str();
    }

    const vector<Argument>& arguments() const { return m_arguments; };
    const int64_t num_argumnets() const { return m_arguments.size(); }
    const Argument& operator[](int index){ return m_arguments[index]; }
};

[[maybe_unused]] std::ostream& operator<<(std::ostream& out, const Command& c){
    out << c.to_string();
    return out;
}

static Command parse_command(){
    Command command;

    cout << "> ";
    std::string user_prompt;
    getline(cin, user_prompt);

    int pos = 0; // position in the user prompt
    while(pos < user_prompt.length() && isspace(user_prompt[pos])) pos++;

    // parse the statement, single word
    int start = pos;
    while(pos < user_prompt.length() && (isalnum(user_prompt[pos]) || user_prompt[pos] == '_')) pos++;
    int end = pos;
//    COUT_DEBUG("start: " << start << ", end: " << end);
    command.m_command = string(user_prompt.c_str() + start, end - start);
    std::transform(command.m_command.begin(), command.m_command.end(), command.m_command.begin(), ::tolower); // all lowercase
//    COUT_DEBUG( "Command: " << command.m_command  );

    while(pos < user_prompt.length()){ // next argument
        if(isspace(user_prompt[pos])){ pos++; continue; } // skip empty spaces at the start

        int string_delimiter = 0;
        if(user_prompt[pos] == '"' || user_prompt[pos] == '\''){
            string_delimiter = user_prompt[pos];
            start = pos +1;
        } else {
            start = pos;
        }
        pos++;

        bool terminated = false;
        while(pos < user_prompt.length() && !terminated){
            if(user_prompt[pos] == '\\'){
                pos+=2;
            } else if (string_delimiter == 0 && isspace(user_prompt[pos])){
                end = pos;
                terminated = true;
            } else if (string_delimiter != 0 && user_prompt[pos] == string_delimiter){
                end = pos; pos++;
                terminated = true;
            } else {
                pos++;
            }
        }
        if(!terminated){ end = pos; }

        if(start < end){
            command.m_arguments.push_back(Argument{user_prompt.c_str() + start, end - start});
        }

        pos ++;
    }

    return command;
}

} // anon

static void run_client_interactive(){
    // Open a connection with the server
    cout << "[client] Connecting to the server at " << cfgclient().get_server_string() << endl;
    network::Client impl { cfgclient().get_server_host(), (int) cfgclient().get_server_port() };

    while(true){
        Command command = parse_command();
        const string& stmt = command.get_command();
        cout << "command parsed: " << command << endl;

        if(stmt == "quit"){
            break;
        } else if(stmt == "terminate"){
            impl.set_terminate_server_on_exit(true);
            break;
        } else if(stmt == "load"){
            if(command.num_argumnets() != 1){
                cout << "[client] [load] invalid number of arguments, expected 1: load (path);" << endl;
            } else {
                impl.load(command[0]);
            }

        } else {
            cout << "[client] ERROR, invalid command: " << stmt << endl;
        }

    }
}

/*****************************************************************************
 *                                                                           *
 *  Main                                                                     *
 *                                                                           *
 *****************************************************************************/

static void run_client(int argc, char* argv[]){
    cout << "[client] Init configuration ... " << endl;
    ClientConfiguration::initialise(argc, argv);

    if(configuration().has_database()){
        cout << "[client] Save the current configuration properties in " << configuration().get_database_path() << endl;
        configuration().save_parameters();
    }

    if(cfgclient().is_interactive()){
        cout << "[client] Interactive mode" << endl;
        run_client_interactive();
    }

    cout << "[client] Done" << endl;
}

int main(int argc, char* argv[]){
    int rc = 0;
    try {
        run_client(argc, argv);
    } catch(common::Error& e){
        cerr << e << endl;
        cerr << "Client terminating due to exception..." << endl;
        rc = 1;
    }

    cfgfree();
    return rc;
}


