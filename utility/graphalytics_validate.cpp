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

#include "graphalytics_validate.hpp"

#include <cctype>
#include <cstdlib>
#include <fstream>
#include <utility>

using namespace std;

namespace utility {

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
#define DEBUG
#define COUT_DEBUG_FORCE(msg) { std::cout << "[GraphalyticsValidate::" << __FUNCTION__ << "] " << msg << std::endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif

/*****************************************************************************
 *                                                                           *
 *  Error                                                                    *
 *                                                                           *
 *****************************************************************************/
#undef CURRENT_ERROR_TYPE
#define CURRENT_ERROR_TYPE ::utility::GraphalyticsValidateError


/*****************************************************************************
 *                                                                           *
 *  Exact match                                                              *
 *                                                                           *
 *****************************************************************************/
void GraphalyticsValidate::exact_match(const std::string& result, const std::string& expected){
    fstream handle_result(result, ios_base::in);
    if(!handle_result.good()) ERROR("The result file does not exist or is not accessible. Path: `"  << result << "'");
    fstream handle_expected(expected, ios_base::in);
    if(!handle_expected.good()) ERROR("The reference file does not exist or is not accessible. Path: `" << expected << "'");

    uint64_t lineno = 0; // current line number

    constexpr size_t buffer_sz = 4096;
    char buffer_result[buffer_sz];
    char buffer_expected[buffer_sz];

    auto parse_distance = [&lineno](char* buffer, const char* buffer_name){
        uint64_t pos = 0;
        char* current = buffer;
        while(pos < buffer_sz && isspace(current[0])) { pos++; current++; };
        if(pos == buffer_sz || current[0] == '\0') ERROR("[lineno=" << lineno << ", file=" << buffer_name << "] The line is empty!");
        if(!isdigit(current[0])) ERROR("[lineno=" << lineno << ", file=" << buffer_name << "] Cannot parse the vertex id in the line `" << buffer << "'");
        char* next = nullptr;
        int64_t vertex_id = strtoll(current, &next, 10);
        current = next;
        while(pos < buffer_sz && isspace(current[0])) { pos++; current++; };
        if(pos == buffer_sz || current[0] == '\0') ERROR("[lineno=" << lineno << ", file=" << buffer_name << "] The line does not contain a distance: `" << buffer << "'");
        if(!isdigit(current[0])) ERROR("[lineno=" << lineno << ", file=" << buffer_name << "] Cannot parse the distance in the line `" << buffer << "'");
        int64_t distance = strtoll(current, nullptr, 10);
        return make_pair(vertex_id, distance);
    };

    while(true){ // read until the eof
        handle_result.getline(buffer_result, buffer_sz);
        handle_expected.getline(buffer_expected, buffer_sz);
        if(handle_result.eof() || handle_expected.eof()) break;

        auto pair_result = parse_distance(buffer_result, "result");
        auto pair_expected = parse_distance(buffer_expected, "expected");

        if(pair_result.first != pair_expected.first){
            ERROR("[lineno=" << lineno << "] VALIDATION ERROR, vertex retrieved: " << pair_result.first << ", vertex expected: " << pair_expected.first);
        } else if (pair_result.second != pair_expected.second){
            ERROR("[lineno=" << lineno << "] VALIDATION ERROR, vertex: " << pair_result.first << " OK, distance retrieved: " << pair_result.first << ", distance expected: " << pair_expected.first);
        }

        lineno++;
    }

    if(handle_expected.eof() && !handle_result.eof()){
        ERROR("[lineno=" << lineno << "] The result file contains more lines [vertices] than then expected/reference output");
    } else if (!handle_expected.eof() && handle_result.eof()){
        ERROR("[lineno=" << lineno << "] The result file contains less lines [vertices] than the expected/reference output");
    }

    handle_result.close();
    handle_expected.close();
}

void GraphalyticsValidate::bfs(const std::string& result, const std::string& expected){
    exact_match(result, expected);
}

} // namespace

