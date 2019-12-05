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

#include "dimacs9_reader.hpp"

#include <cassert>
#include <cctype>
#include "graph/edge.hpp"
#include "configuration.hpp"
#include "utility.hpp"

using namespace std;

#undef CURRENT_ERROR_TYPE
#define CURRENT_ERROR_TYPE ::gfe::reader::ReaderError

namespace gfe::reader {

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG
#define COUT_DEBUG_FORCE(msg) { std::cout << "[Dimacs9Reader::" << __FUNCTION__ << "] " << msg << std::endl; }
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
Dimacs9Reader::Dimacs9Reader(const std::string& path) : m_handle(init_fstream(path)) {
    m_buffer[sizeof(m_buffer) -1] = '\0'; // avoid overflows
}

Dimacs9Reader::~Dimacs9Reader(){
    m_handle.close();
}

bool Dimacs9Reader::is_directed() const {
    return true;
}

/*****************************************************************************
 *                                                                           *
 *  Parser                                                                   *
 *                                                                           *
 *****************************************************************************/
bool Dimacs9Reader::parse_line(char*& line){
    if(line == nullptr) return false;
    if(line[0] != 'a'){ return false; } // this does not represent an edge
    do{ line++; } while(isspace(line[0])); // ignore the following spaces
    return true;
}

bool Dimacs9Reader::is_number(const char* marker){
    return marker != nullptr && (marker[0] >= '0' && marker[0] <= '9');
}

bool Dimacs9Reader::read(graph::WeightedEdge& edge) {
    if(!m_handle.good()) return false;

    // read the next line that is not a comment
    bool skip { true } ;
    char* current = nullptr;
    do {
        m_handle.getline(m_buffer, sizeof(m_buffer) -1);
        current = m_buffer;
        skip = ! parse_line(current);
#if defined(DEBUG)
      if(skip) { COUT_DEBUG("line `" << m_buffer << "' skipped"); }
#endif
    } while (skip && m_handle.good());
    if(skip) return false; // ended with a comment

    // read the source vertex
    char* next { nullptr };
    if(!is_number(current)) ERROR("line: `" << m_buffer << "', cannot read the source vertex");
    edge.m_source = strtoull(current, &next, /* base */ 10);

    // read the destination vertex
    while(isspace(next[0])) next++;
    current = next;
    if(!is_number(current)) ERROR("line: `" << m_buffer << "', cannot read the destination vertex");
    edge.m_destination = strtoull(current, &next, 10);

    // read the edge weight
    while(isspace(next[0])) next++;
    current = next;
    if(!is_number(current)) ERROR("line: `" << m_buffer << "', cannot read the weight");
    edge.m_weight = strtoul(current, nullptr, 10);

    COUT_DEBUG("edge parsed: " << edge);

    return true;
}

} // namespace


