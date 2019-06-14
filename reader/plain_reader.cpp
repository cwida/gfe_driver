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

#include "plain_reader.hpp"

#include <fstream>
#include <string>
#include "graph/edge.hpp"
#include "configuration.hpp"
#include "utility.hpp"

using namespace std;

#undef CURRENT_ERROR_TYPE
#define CURRENT_ERROR_TYPE reader::ReaderError

namespace reader {

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG
#define COUT_DEBUG_FORCE(msg) { std::cout << "[PlainReader::" << __FUNCTION__ << "] " << msg << std::endl; }
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
PlainReader::PlainReader(const std::string& path, bool is_weighted) : PlainReader(path, is_weighted, configuration().max_weight()) { }

PlainReader::PlainReader(const string& path, bool is_weighted, uint32_t max_weight) :
        m_handle(init_fstream(path)), m_is_weighted(is_weighted), m_max_weight(max_weight),
        m_random_generator(configuration().seed() + 130233320) {
    LOG("[PlainWeightedReader] Reading `" << path << "' ...");
    COUT_DEBUG("path: " << path << ", is_weighted: " << m_is_weighted << ", max_weight: " << m_max_weight);
}

PlainReader::~PlainReader(){
    m_handle.close();
}

/*****************************************************************************
 *                                                                           *
 *  Parser                                                                   *
 *                                                                           *
 *****************************************************************************/

bool PlainReader::ignore_line(const char* line){
    if(line == nullptr) return true;
    while(isspace(line[0])) line++;
    return line[0] == '#' || line[0] == '\0';
}

bool PlainReader::ignore_line(const std::string& line){
    return ignore_line(line.c_str());
}

bool PlainReader::is_number(const char* marker){
    return marker != nullptr && (marker[0] >= '0' && marker[0] <= '9');
}

bool PlainReader::read(graph::WeightedEdge& edge) {
    if(!m_handle.good()) return false;

    // read the next line that is not a comment
    string current_line;
    bool skip { true } ;
    do {
        getline(m_handle, current_line);

        skip = ignore_line(current_line);
#if defined(DEBUG)
      if(skip) { COUT_DEBUG("line: `" << current_line << "' is a comment or an empty line, skipped"); }
#endif
    } while (skip && m_handle.good());
    if(skip) return false; // ended with a comment

    // read the source
    char* next { nullptr };
    const char* current = current_line.c_str();
    while(isspace(current[0])) current++;
    if(!is_number(current)) ERROR("line: `" << current_line << "', cannot read the source vertex");
    edge.m_source = strtoull(current, &next, /* base */ 10);

    while(isspace(next[0])) next++;
    current = next;
    if(!is_number(current)) ERROR("line: `" << current_line << "', cannot read the destination vertex");
    edge.m_destination = strtoull(current, &next, 10);

    if(m_is_weighted){
        while(isspace(next[0])) next++;
        current = next;
        if(!is_number(current)) ERROR("line: `" << current_line << "', cannot read the weight");
        edge.m_weight = strtoul(current, nullptr, 10);
    } else if (m_max_weight == 1){
        edge.m_weight = 1;
    } else { // provide a random weight
        uniform_int_distribution<uint32_t> distribution{1, static_cast<uint32_t>(m_max_weight)};
        edge.m_weight = distribution(m_random_generator);
    }
    COUT_DEBUG("edge parsed: " << edge);

    return true;
}

} // namespace reader
