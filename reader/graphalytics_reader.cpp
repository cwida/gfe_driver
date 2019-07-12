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

#include "graphalytics_reader.hpp"

#include <cassert>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <regex>
#include "common/filesystem.hpp"
#include "graph/edge.hpp"
#include "configuration.hpp"
#include "utility.hpp"

using namespace std;

#undef CURRENT_ERROR_TYPE
#define CURRENT_ERROR_TYPE ::reader::ReaderError

namespace reader {

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG
#define COUT_DEBUG_FORCE(msg) { std::cout << "[GraphalyticsReader::" << __FUNCTION__ << "] " << msg << std::endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif

/*****************************************************************************
 *                                                                           *
 *  Initialisation                                                           *
 *                                                                           *
 *****************************************************************************/
GraphalyticsReader::GraphalyticsReader(const std::string& path_properties){
    if(!common::filesystem::file_exists(path_properties)) ERROR("THe given file does not exist: " << path_properties);
    string abs_path_properties = common::filesystem::absolute_path(path_properties);
    m_properties.insert({string("property-file"), abs_path_properties});
    string basedir = common::filesystem::directory(abs_path_properties);

    regex pattern { "^\\s*graph\\.([A-Za-z0-9_.-]+?)\\s*=\\s*([^#\n]+?)\\s*" };

    string name;

    COUT_DEBUG("Parsing the property file: " << path_properties);
    fstream handle { path_properties.c_str() };
    while(handle.good()){
        string line;
        getline(handle, line);

        smatch matches;
        if ( regex_match(line, matches, pattern) ) { // side effect => populate matches
            string key = matches[1];
            string value = matches[2];

            if(name.empty()){
                size_t pos = key.find('.');
                if(pos == string::npos) ERROR("Cannot parse the name of the graph (expected to find graph.<GRAPH_NAME>.<PROPERTY_NAME>), key=graph." << key);
                name = key.substr(0, pos);
                COUT_DEBUG("name=" << name);
                m_properties.insert({string("name"), name});
            }

            if( key.substr(0, name.length()) != name ){
                LOG("[GraphalyticsReader] Warning, line skipped `" << line << "', graph name does not match what expected: " << name);
                continue;
            }
            key = key.substr(name.length() +1);

            if(key == "vertex-file" || key == "edge-file") {
                if(value.empty()) ERROR("Empty path specified for the key " << key << ", line: " << line);
                if(value[0] != '/'){
                    value = common::filesystem::absolute_path(basedir + "/" + value);
                }
            } else if (key == "directed"){
                regex pattern_yes { "^\\s*(yes|true|1)\\s*$", regex_constants::icase };
                regex pattern_no { "^\\s*(no|false|0)\\s*$", regex_constants::icase };

                if(regex_search(value, pattern_yes)){
                    COUT_DEBUG("The graph is directed");
                    m_directed = true;
                } else if (regex_search(value, pattern_no)){
                    COUT_DEBUG("The graph is undirected");
                    m_directed = false;
                } else {
                    ERROR("Cannot determine whether the graph is directed or not. The property value is `" << value << "'");
                }
            }

            COUT_DEBUG("key: " << key << ", value: " << value);
            m_properties.insert({key, value});
        }
    }
    handle.close();

    // check that the key vertex-file, edge-file and directed are present
    if(m_properties.find("vertex-file") == m_properties.end()) ERROR("The property `vertex-file' is not set in the property file");
    if(m_properties.find("edge-file") == m_properties.end()) ERROR("The property `edge-file' is not set in the property file");
    COUT_DEBUG("vertex-file: " << get_path_vertex_list() << ", edge-file: " << get_path_edge_list() << ", is_directed: " << is_directed());
}

GraphalyticsReader::~GraphalyticsReader(){
    close();
}


/*****************************************************************************
 *                                                                           *
 *  Properties                                                               *
 *                                                                           *
 *****************************************************************************/

string GraphalyticsReader::get_property(const std::string& key) const {
    auto item = m_properties.find(key);
    return (item == m_properties.end()) ? "" : item->second;
}

string GraphalyticsReader::get_path_vertex_list() const {
    auto item = m_properties.find("vertex-file");
    if(item == m_properties.end()) ERROR("The property `vertex-file' is not set");
    return item->second;
}

string GraphalyticsReader::get_path_edge_list() const {
    auto item = m_properties.find("edge-file");
    if(item == m_properties.end()) ERROR("The property `edge-file' is not set");
    return item->second;
}

bool GraphalyticsReader::is_directed() const {
    return m_directed;
}

bool GraphalyticsReader::is_weighted() const {
    return m_is_weighted;
}

void GraphalyticsReader::set_emit_directed_edges(bool value){
    m_emit_directed_edges = value;
}

/*****************************************************************************
 *                                                                           *
 *  Handles                                                                  *
 *                                                                           *
 *****************************************************************************/

static void handle_close(void*& ptr_handle){
    if(ptr_handle == nullptr) return; // nop
    auto handle = reinterpret_cast<fstream*>(ptr_handle);
    handle->close();
    delete handle;
    ptr_handle = nullptr;
}

void GraphalyticsReader::close(){
    handle_close(m_handle_edge_file);
    handle_close(m_handle_vertex_file);

}

void GraphalyticsReader::reset(){
    // read_edge and read_vertex are going to reinit the handle if they are called again
    close();
    m_last_reported = true;
}

/*****************************************************************************
 *                                                                           *
 *  Read                                                                     *
 *                                                                           *
 *****************************************************************************/

bool GraphalyticsReader::read(graph::WeightedEdge& edge){
    return read_edge(edge);
}

bool GraphalyticsReader::read_edge(graph::WeightedEdge& edge){
    if(m_handle_edge_file == nullptr) {
        COUT_DEBUG("Opening the input stream `" << get_path_edge_list() << "'");
        m_handle_edge_file = new fstream(get_path_edge_list());
    }
    fstream* handle = reinterpret_cast<fstream*>(m_handle_edge_file);
    assert(handle != nullptr && "Null pointer");

    if(!is_directed() && m_emit_directed_edges && !m_last_reported){
        std::swap(m_last_source, m_last_destination);
        m_last_reported = true;
    } else {
        if(!handle->good()) return false;

        // read the next line that is not a comment
        string current_line;
        bool skip { true } ;
        do {
            getline(*handle, current_line);

            skip = ignore_line(current_line);
#if defined(DEBUG)
          if(skip) { COUT_DEBUG("line: `" << current_line << "' is a comment or an empty line, skipped"); }
#endif
        } while (skip && handle->good());
        if(skip) return false; // ended with a comment
        COUT_DEBUG("Parse line: `" << current_line << "'");

        // read the source
        char* next { nullptr };
        const char* current = current_line.c_str();
        while(isspace(current[0])) current++;
        if(!is_number(current)) ERROR("line: `" << current_line << "', cannot read the source vertex");
        m_last_source = strtoull(current, &next, /* base */ 10);

        while(isspace(next[0])) next++;
        current = next;
        if(!is_number(current)) ERROR("line: `" << current_line << "', cannot read the destination vertex");
        m_last_destination = strtoull(current, &next, 10);

        if(is_weighted()){
            while(isspace(next[0])) next++;
            current = next;
            if(!is_number(current)) ERROR("line: `" << current_line << "', cannot read the weight");
            m_last_weight = strtod(current, nullptr);
        } else {
            m_last_weight = 1.0;
        }

        m_last_reported = false;
    }

    edge.m_source = m_last_source;
    edge.m_destination = m_last_destination;
    edge.m_weight = m_last_weight;
    COUT_DEBUG("edge parsed: " << edge);

    return true;
}

bool GraphalyticsReader::read_vertex(uint64_t& out_vertex){
    out_vertex = 0; // init

    if(m_handle_vertex_file == nullptr) {
        COUT_DEBUG("Opening the input stream `" << get_path_vertex_list() << "'");
        m_handle_vertex_file = new fstream(get_path_vertex_list());
    }
    fstream* handle = reinterpret_cast<fstream*>(m_handle_vertex_file);
    assert(handle != nullptr && "Null pointer");

    if(!handle->good()) return false;

    // read the next line that is not a comment
    string current_line;
    bool skip { true } ;
    do {
        getline(*handle, current_line);

        skip = ignore_line(current_line);
#if defined(DEBUG)
      if(skip) { COUT_DEBUG("line: `" << current_line << "' is a comment or an empty line, skipped"); }
#endif
    } while (skip && handle->good());
    if(skip) return false; // ended with a comment
    COUT_DEBUG("Parse line: `" << current_line << "'");

    const char* current = current_line.c_str();
    while(isspace(current[0])) current++;
    if(!is_number(current)) ERROR("line: `" << current_line << "', cannot read the vertex id");
    out_vertex = strtoull(current, nullptr, /* base */ 10);

    return true;
}


bool GraphalyticsReader::ignore_line(const char* line){
    if(line == nullptr) return true;
    while(isspace(line[0])) line++;
    return line[0] == '#' || line[0] == '\0';
}

bool GraphalyticsReader::ignore_line(const std::string& line){
    return ignore_line(line.c_str());
}

bool GraphalyticsReader::is_number(const char* marker){
    return marker != nullptr && (marker[0] >= '0' && marker[0] <= '9');
}


} // namespace

