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

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "common/error.hpp"

using namespace std;

namespace gfe::utility {

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG
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
#define CURRENT_ERROR_TYPE ::gfe::utility::GraphalyticsValidateError
#define FATAL ERROR
#define ERROR_INIT uint64_t error_count = 0;
#define ERROR_COUNT(msg) if(max_num_errors == 1){ FATAL(msg); } else { \
	error_count ++; /* error_count is a local var that should be already defined, init to 0 */ \
	std::cerr << "VALIDATION ERROR #" << error_count << ": " << msg << endl; \
	if(max_num_errors > 1 && error_count == max_num_errors) { \
		FATAL("Reached the maximum number of errors: " << max_num_errors); \
	} \
}
#define ERROR_EXIT if(error_count > 0){ FATAL("Validation found " << error_count << " mismatches"); }

/*****************************************************************************
 *                                                                           *
 *  Helpers                                                                  *
 *                                                                           *
 *****************************************************************************/
constexpr size_t BUFFER_SZ = 4096;

template<typename T> static bool validate_value_typed(char* buffer){ return false; /* it should never be instantiated */ }
template<> bool validate_value_typed<int64_t>(char* buffer){ return isdigit(buffer[0]); }
template<> bool validate_value_typed<double>(char* buffer){ return isdigit(buffer[0]) || buffer[0] == '.' || buffer[0] == 'i' /* infinity */ ; }
template<typename T> static T parse_value_typed(char* buffer){ }
template<> int64_t parse_value_typed<int64_t>(char* buffer){ return strtoll(buffer, nullptr, 10); }
template<> double parse_value_typed<double>(char* buffer){ return strtod(buffer, nullptr); }

template<typename T>
static pair<int64_t, T> parse_value(uint64_t lineno, char* buffer, const char* buffer_name){
    uint64_t pos = 0;
    char* current = buffer;
    while(pos < BUFFER_SZ && isspace(current[0])) { pos++; current++; };
    if(pos == BUFFER_SZ || current[0] == '\0') FATAL("[lineno=" << lineno << ", file=" << buffer_name << "] The line is empty!");
    if(!isdigit(current[0])) FATAL("[lineno=" << lineno << ", file=" << buffer_name << "] Cannot parse the vertex id in the line `" << buffer << "'");
    char* next = nullptr;
    int64_t vertex_id = strtoll(current, &next, 10);
    current = next;
    while(pos < BUFFER_SZ && isspace(current[0])) { pos++; current++; };
    if(pos == BUFFER_SZ || current[0] == '\0') FATAL("[lineno=" << lineno << ", file=" << buffer_name << "] The line does not contain a value: `" << buffer << "'");
    if(!validate_value_typed<T>(current)) FATAL("[lineno=" << lineno << ", file=" << buffer_name << "] Cannot parse the value in the line `" << buffer << "'");
    T value = parse_value_typed<T>(current);
    return make_pair(vertex_id, value);
};

/**
 * Relabel the vertex ID according to the given map
 */
static pair<int64_t, int64_t> relabel(const pair<int64_t, int64_t>& tuple, uint64_t lineno, const std::unordered_map<uint64_t, uint64_t>* vtx_map, bool relabel_value) {
    if(vtx_map == nullptr) return tuple; // nothing to relabel

    pair<int64_t, int64_t> result = tuple;

    { // restrict the scope
        auto remap = vtx_map->find(tuple.first); // vertex ID
        if(remap == vtx_map->end()){
            FATAL("[lineno=" << lineno << "] VALIDATION ERROR, cannot remap the vertex ID `" << tuple.first << "'");
        }
        result.first = remap->second;
    }

    if(relabel_value){
        auto remap = vtx_map->find(tuple.second); // value
        if(remap == vtx_map->end()){
            FATAL("[lineno=" << lineno << "] VALIDATION ERROR, cannot remap the value for the vertex `" << tuple.second << "'");
        }
        result.second = remap->second;
    }

    return result;
}

namespace { template<typename T> struct Tuple { int64_t vertex_id; T value; uint64_t lineno; }; }

/**
 * Read the content of the given reference/expected file, return an hash map vertex id -> <value, line in the file>
 */
template<typename T>
static unordered_map</* vertex id */ int64_t, /* value */ Tuple<T>> read_results(const std::string& path_to_file){
    fstream handle(path_to_file, ios_base::in);
    if(!handle.good()) FATAL("The result file does not exist or is not accessible. Path: `"  << path_to_file << "'");

    uint64_t lineno = 0; // current line number
    char buffer[BUFFER_SZ]; // read the current line from the file
    unordered_map<int64_t, Tuple<T>> result;

    lineno = 0;
    while(true){
        handle.getline(buffer, BUFFER_SZ);
        if(handle.eof()) break;
        auto v = parse_value<T>(lineno, buffer, "result");
        auto rc = result.insert({ v.first, Tuple<T>{v.first, v.second, lineno} });
        if(!rc.second){ FATAL("[lineno=" << lineno << ", file=" << path_to_file << "] The vertex " << v.first << " is a duplicate, already defined at line #" << rc.first->second.lineno); }
        lineno++;
    }

    handle.close();

    return result;
}

/*****************************************************************************
 *                                                                           *
 *  Exact match                                                              *
 *                                                                           *
 *****************************************************************************/
void GraphalyticsValidate::exact_match(const std::string& path_result, const std::string& path_expected, uint64_t max_num_errors, const vertex_map_t* vtx_map, bool vtx_relabel_values){
    ERROR_INIT

    fstream handle_expected(path_expected, ios_base::in);
    if(!handle_expected.good()) FATAL("The reference file does not exist or is not accessible. Path: `" << path_expected << "'");
    auto hash_results = read_results<int64_t>(path_result); // hash map: vertex id -> <value, line>

    // now parse `reference'
    char buffer[BUFFER_SZ];
    uint64_t lineno = 0;
    while(true){
        handle_expected.getline(buffer, BUFFER_SZ);
        if(handle_expected.eof()) break;
        auto t_expected = relabel( parse_value<int64_t>(lineno, buffer, "expected"), lineno, vtx_map, vtx_relabel_values );

        if(lineno >= hash_results.size()){
            ERROR_COUNT("[lineno=" << lineno << "] VALIDATION ERROR, the reference contains more vertices than the actual result file");
        } else {
            const auto& t_result_it = hash_results.find(t_expected.first);
            if (t_result_it == end(hash_results)){
                ERROR_COUNT("[line number reference: " << lineno << "] VALIDATION ERROR, the vertex " << t_expected.first << " is present in the reference (" << path_expected << ") but not in the results (" << path_result << ") ");
            } else if (t_result_it->second.vertex_id != t_expected.first){ // this should never occur, as the index of the hash map is the vertex id itself
                ERROR_COUNT("[line number result: " << t_result_it->second.lineno << ", reference: " << lineno << "] VALIDATION ERROR, vertex retrieved: " << t_result_it->second.vertex_id << ", vertex expected: " << t_expected.first << " ");
            } else if (t_expected.second != t_result_it->second.value){
                ERROR_COUNT("[line number result: " << t_result_it->second.lineno << ", reference: " << lineno << "] VALIDATION ERROR, vertex: " << t_result_it->second.vertex_id << " matches, but value retrieved: " << t_result_it->second.value << " != value expected: " << t_expected.second);
            }
        }

        lineno++;
    }

    if(lineno < hash_results.size()){
    	ERROR_COUNT("The result file contains more lines [vertices] than the expected/reference output. Vertices in the result file: " << hash_results.size() << ", vertices expected: " << lineno);
    }

    handle_expected.close();

    ERROR_EXIT
}

void GraphalyticsValidate::bfs(const std::string& result, const std::string& expected, uint64_t max_num_errors, const vertex_map_t* vertex_map){
    exact_match(result, expected, max_num_errors, vertex_map, false);
}

void GraphalyticsValidate::cdlp(const std::string& result, const std::string& expected, uint64_t max_num_errors, const vertex_map_t* vertex_map){
    exact_match(result, expected, max_num_errors, vertex_map, true);
}

/*****************************************************************************
 *                                                                           *
 *  Epsilon match                                                            *
 *                                                                           *
 *****************************************************************************/
void GraphalyticsValidate::epsilon_match(const std::string& path_result, const std::string& path_expected, double epsilon, uint64_t max_num_errors, const vertex_map_t* vertex_map){
    ERROR_INIT

    fstream handle_expected(path_expected, ios_base::in);
    if(!handle_expected.good()) FATAL("The reference file does not exist or is not accessible. Path: `" << path_expected << "'");
    auto hash_results = read_results<double>(path_result); // hash map: vertex id -> <value, line>

    uint64_t lineno = 0; // current line number
    char buffer[BUFFER_SZ]; // current line read from the file

    while(true){ // read until the eof
        handle_expected.getline(buffer, BUFFER_SZ);
        if(handle_expected.eof()) break;

        auto t_expected = parse_value<double>(lineno, buffer, "expected");
        if(vertex_map != nullptr) { t_expected.first = vertex_map->find(t_expected.first)->second; } // remap the vertex

        if(lineno >= hash_results.size()){
            ERROR_COUNT("[lineno=" << lineno << "] VALIDATION ERROR, the reference contains more vertices than the actual result file");
        } else {

            const auto& t_result_it = hash_results.find(t_expected.first);
            if (t_result_it == end(hash_results)){
                ERROR_COUNT("[line number reference: " << lineno << "] VALIDATION ERROR, the vertex " << t_expected.first << " is present in the reference (" << path_expected << ") but not in the results (" << path_result << ") ");
            } else if (t_result_it->second.vertex_id != t_expected.first){ // this should never occur, as the index of the hash map is the vertex id itself
                ERROR_COUNT("[line number result: " << t_result_it->second.lineno << ", reference: " << lineno << "] VALIDATION ERROR, vertex retrieved: " << t_result_it->second.vertex_id << ", vertex expected: " << t_expected.first << " ");
            } else {
                double value_result = t_result_it->second.value;
                double value_expected = t_expected.second;


                double error = abs(value_result - value_expected) / value_expected;
                COUT_DEBUG("vertex: " << t_result_it->second.vertex_id << ", value: " << value_result << ", expected: " << value_expected << ", error: " << error);
                if (error > epsilon){
                    ERROR_COUNT("[lineno result: " << t_result_it->second.lineno << ", reference:" << lineno << "] VALIDATION ERROR, vertex: " << t_result_it->second.vertex_id << " matches, but "
                            "value retrieved: " << value_result << ", value expected: " << value_expected << ", error: " << error << ", tolerance (epsilon): " << epsilon);
                }
            }
        }

        lineno++;
    }

    if(lineno < hash_results.size()){
        ERROR_COUNT("The result file contains more lines [vertices] than the expected/reference output. Vertices in the result file: " << hash_results.size() << ", vertices expected: " << lineno);
    }

    handle_expected.close();

    ERROR_EXIT
}

void GraphalyticsValidate::pagerank(const std::string& result, const std::string& expected, uint64_t max_num_errors, const vertex_map_t* vertex_map){
    epsilon_match(result, expected, /* as defined in LDBC Graphalytics ~ PageRankValidationTest.java */ 0.0001, max_num_errors, vertex_map);
}

void GraphalyticsValidate::lcc(const std::string& result, const std::string& expected, uint64_t max_num_errors, const vertex_map_t* vertex_map){
    epsilon_match(result, expected, 0.0001, max_num_errors, vertex_map);
//    epsilon_match(result, expected, /* as defined in LDBC Graphalytics ~ LocalClusteringCoefficientValidationTest.java */ 0.000001);
}

void GraphalyticsValidate::sssp(const std::string& result, const std::string& expected, uint64_t max_num_errors, const vertex_map_t* vertex_map){
    epsilon_match(result, expected, 0.0001, max_num_errors, vertex_map);
}


/*****************************************************************************
 *                                                                           *
 *  Equivalence match                                                        *
 *                                                                           *
 *****************************************************************************/

void GraphalyticsValidate::equivalence_match(const std::string& result, const std::string& expected, uint64_t max_num_errors, const vertex_map_t* vertex_map){
    ERROR_INIT

    fstream handle_result(result, ios_base::in);
    if(!handle_result.good()) FATAL("The result file does not exist or is not accessible. Path: `"  << result << "'");
    fstream handle_expected(expected, ios_base::in);
    if(!handle_expected.good()) FATAL("The reference file does not exist or is not accessible. Path: `" << expected << "'");

    char buffer[BUFFER_SZ]; // current line being processed
    uint64_t lineno = 0; // current line number
    unordered_map<int64_t, int64_t> components_result; // read the results
    unordered_map<int64_t, int64_t> ref_mapping; // component[ref] -> component[exp]
    unordered_set<int64_t> unique_components; // check that component[exp] does not belong to different components in ref

    // process the file with the results
    while(true){
        handle_result.getline(buffer, BUFFER_SZ);
        if(handle_result.eof()) break;
        auto t = parse_value<int64_t>(lineno, buffer, "result");
        components_result[t.first] = t.second;
        lineno++;
    }

    handle_result.close();

    // process the reference file
    lineno = 0;
    while(true){
        handle_expected.getline(buffer, BUFFER_SZ);
        if(handle_expected.eof()) break;
        auto t_ref = parse_value<int64_t>(lineno, buffer, "reference");
        if(vertex_map != nullptr) { t_ref.first = vertex_map->find(t_ref.first)->second; } // relabel the vertex

        // first of all, does this vertex exist in the result file?
        auto ptr_t_res = components_result.find(t_ref.first);
        if(ptr_t_res == end(components_result)){
            ERROR_COUNT("[lineno reference:" << lineno << "] VALIDATION ERROR, the vertex " << t_ref.first << " is expected but not present in the result file");
        } else {
            auto t_res = *ptr_t_res;

            // is the first time we see this vertex in the ref file?
            auto mapping = ref_mapping.insert({t_ref.second, t_res.second});
            if(mapping.second){ // as we have seen for the first time this component in ref, check we haven't seen previously also in exp
                auto insert_rc = unique_components.insert(t_res.second);
                if(!insert_rc.second){ // we were not able to insert `result' in the list of unique_values => a mapping for `result' already exists
                    ERROR_COUNT("[lineno reference:" << lineno << "] VALIDATION ERROR, vertex: " << t_ref.first << ", the component " << t_res.second << " is associated to a single component in the result file but "
                            "belongs to two different components in the reference file");
                }
            } else if(!mapping.second && mapping.first->second != t_res.second) { // this mapping already exists, but the two components don't match
                ERROR_COUNT("[lineno reference:" << lineno << "] VALIDATION ERROR, vertex: " << t_ref.first << ", invalid mapping, component in the result file: " << t_res.second <<
                        ", expected value: " << mapping.first->second << " (in ref. file, mapped to value: " << t_ref.second << ")");
            }

        }

        lineno++;
    }

    if(lineno < components_result.size()){
        ERROR_COUNT("The result file contains more lines [vertices] than the expected/reference output. Vertices result:  " << components_result.size() << ", vertices expected: " << lineno);
    }

    handle_expected.close();

    ERROR_EXIT
}

void GraphalyticsValidate::wcc(const std::string& result, const std::string& expected, uint64_t max_num_errors, const vertex_map_t* vertex_map){
    equivalence_match(result, expected, max_num_errors, vertex_map);
}

} // namespace
