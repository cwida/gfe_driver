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

#include "metis_reader.hpp"

#include <cassert>
#include <cctype>
#include <regex>
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
#define COUT_DEBUG_FORCE(msg) { std::cout << "[MetisReader::" << __FUNCTION__ << "] " << msg << std::endl; }
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

MetisReader::MetisReader(const std::string& path) : m_handle(init_fstream(path)), m_random_generator(configuration().seed() + 12908478){
    LOG("[MetisReader] Reading `" << path << "' ...");

    // parse the header
    fetch_next_line();
    if(m_handle.bad()) ERROR("Cannot read the header");
    regex pattern { "^\\s*(\\d+)\\s+(\\d+)(?:\\s+([01]{3})(?:\\s+(\\d+))?)?\\s*$" };
    smatch matches;
    if ( ! regex_match(m_current_line, matches, pattern) ) { // side effect => populate matches
        ERROR("[MetisReader] Cannot parse the header: `" << m_current_line << "'");
    }
    assert(matches.ready() && "Matches should have been populated by ::regex_match");
    COUT_DEBUG("Number of matches: " << matches.size());
    // matches[0] is the entire string matched
    m_num_vertices = stoull(matches[1]);
//    m_num_edges = stoull(matches[2]); // ignore this field, we're not going to validate the number of edges parsed

    if(!matches[3].matched){ // optional field, according to the manual, it means that all weights not explicitly provided but are implicitly set to 1
        m_num_vertex_weights = 0;
        m_has_edge_weight = false;
        m_is_edge_weight_constant = true;
        m_edge_weight = 1;
    } else { // the format is explicitly present. It consists of three digits in {0, 1}
        string format = matches[3];
        assert(format.length() == 3 && R"(it should consist of a string \d\d\d, no one digit less, no one digit more)");

        if(matches[4].matched){ // number of vertex weights to skip
            int ncon = stoi(matches[4]);
            assert((ncon != 0 || format.c_str()[1] == '0') && "ncon == 0 => vertex_weights == 0"); // !A or B
            assert((ncon == 0 || format.c_str()[1] == '1') && "ncon != 0 => vertex_weights == 1"); // !A or B
            m_num_vertex_weights = /* vertex size */ (format.c_str()[0] == '1') + ncon;
        } else {
            m_num_vertex_weights = /* vertex size */ (format.c_str()[0] == '1') + (format.c_str()[1] == '1');
        }

        m_has_edge_weight = format.c_str()[2] == '1';
    }
    // reset the ptr used by #read_value
    m_current_line_read_value_ptr = nullptr;

    COUT_DEBUG("has edge weights: " << m_has_edge_weight << ", num vertex weights: " << m_num_vertex_weights);
}

MetisReader::~MetisReader(){
    m_handle.close();
}

bool MetisReader::is_directed() const {
    return true;
}

/*****************************************************************************
 *                                                                           *
 *  Parser                                                                   *
 *                                                                           *
 *****************************************************************************/

bool MetisReader::fetch_next_line(){
    if(!m_handle.good()) return false;
    bool comment = false;
    do {
        getline(m_handle, m_current_line);
        m_lineno++;

       comment = is_comment();
#if defined(DEBUG)
      if(comment) { COUT_DEBUG("[" << m_lineno << "] line: `" << m_current_line << "' is a comment, skipped"); }
#endif
    } while (comment && m_handle.good());

    if(comment) return false; // ended with a comment

    m_current_line_read_value_ptr = (char*) m_current_line.c_str();
    while(isspace(m_current_line_read_value_ptr[0])) m_current_line_read_value_ptr++;

    return true;
}

bool MetisReader::is_comment() const {
    static regex pattern { "^\\s*%" };
    return regex_search(m_current_line, pattern);
}

std::pair<bool, uint64_t> MetisReader::read_value(){
    std::pair<bool, uint64_t> result { /* were we able to read a value ? */ false, /* the value read */ 0 };

    if(m_current_line_read_value_ptr != nullptr && (m_current_line_read_value_ptr[0] >= '0' && m_current_line_read_value_ptr[0] <= '9')){
        char* next = nullptr;
        result.first = true;
        result.second = strtoull(m_current_line_read_value_ptr, &next, 10);
        while(isspace(next[0])) next++;
        m_current_line_read_value_ptr = next;
    }

    return result;
}

bool MetisReader::read(graph::WeightedEdge& e) {
    // read the next destination vertex
    auto dest_vertex = read_value();

    if(!dest_vertex.first) { // move to the next line that is not a comment
        if (! fetch_next_line() ) return false; // side effect, it updates the value of m_current_line
        COUT_DEBUG("Parse next line: " << m_current_line);

        if(m_handle.bad()) return false;
        m_edge_vertex1++; // the first vertex starts from 1 (not zero!)
        if(m_edge_vertex1 > m_num_vertices) {
            assert(m_edge_vertex1 > 0 && "m_edge_vertex -1 will cause an underflow");
            ERROR("line no: " << m_lineno << ", line: `" << m_current_line << "', constraint not respected: the number of vertices in the header ("
                    << m_num_vertices << ") is less than the number of vertices read: " << (m_edge_vertex1 -1));
        }

        // skip the vertex weights
        for(int i = 0; i < m_num_vertex_weights; i++){
            auto result = read_value();
            if(!result.first) {
                ERROR("line no: " << m_lineno << ", vertex: " << m_edge_vertex1 << ", cannot skip the first " << m_num_vertex_weights << " vertex weight. Stuck at: " << i << "/" << m_num_vertex_weights);
            }
        }

        dest_vertex = read_value();
        if(!dest_vertex.first) return false; // we're done
    }

    COUT_DEBUG("Current line: " << m_current_line);

    m_edge_vertex2 = dest_vertex.second;

    // edge weight
    if(m_is_edge_weight_constant){
        assert(m_has_edge_weight == false && "The edge weight should not be explicitly provided");
        // nop
    } else if(m_has_edge_weight){
        auto info_weight = read_value();
        if(!info_weight.first) {
            ERROR("line no: " << m_lineno << ", missing weight for the edge " << m_edge_vertex1 << " -> " << m_edge_vertex2 << ", line: " << m_current_line);
        } else if(info_weight.second > std::numeric_limits<int32_t>::max()) {
            ERROR("line no: " << m_lineno << ", the weight for the edge " << m_edge_vertex1 << " -> " << m_edge_vertex2 << " "
                    "is too big: " << info_weight.second << ", line: " << m_current_line);
        }
        m_edge_weight = info_weight.second;
    } else { // assign an edge weight randomly
        uniform_int_distribution<uint32_t> distribution{1, static_cast<uint32_t>(configuration().max_weight())};
        m_edge_weight = distribution(m_random_generator);
    }

    COUT_DEBUG("parsed edge: " << m_edge_vertex1 << " -> " << m_edge_vertex2 << ", weight: " << m_edge_weight);

    e.m_source = m_edge_vertex1;
    e.m_destination = m_edge_vertex2;
    e.m_weight = m_edge_weight;

    return true;
}

} // namespace

