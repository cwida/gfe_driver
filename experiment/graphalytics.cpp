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

#include "graphalytics.hpp"

#include <algorithm>
#include <cstdio> // mkdtemp
#include <cstring>
#include <filesystem>
#include <iostream>
#include <string>
#include <sstream>

#include "common/database.hpp"
#include "common/filesystem.hpp"
#include "common/timer.hpp"
#include "library/interface.hpp"
#include "reader/graphalytics_reader.hpp"
#include "utility/graphalytics_validate.hpp"
#include "configuration.hpp"
#include "statistics.hpp"

using namespace common;
using namespace gfe::utility;
using namespace std;

namespace gfe::experiment {

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
#define DEBUG
#define COUT_DEBUG_FORCE(msg) { std::cout << "[Graphalytics @ " << __FUNCTION__ << "] " << msg << std::endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif


/*****************************************************************************
 *                                                                           *
 *  GraphalyticsAlgorithms                                                   *
 *                                                                           *
 ****************************************************************************/
GraphalyticsAlgorithms::GraphalyticsAlgorithms(const string& path) {
    COUT_DEBUG("Read the properties from path: " << path);
    reader::GraphalyticsReader props { path };

    stringstream ss { props.get_property("algorithms") };
    string algorithm_name;
    while (getline(ss, algorithm_name, ',')){
        algorithm_name.erase(std::remove_if(begin(algorithm_name), end(algorithm_name), ::isspace), end(algorithm_name));
        transform(algorithm_name.begin(), algorithm_name.end(), algorithm_name.begin(), ::tolower);
        if(algorithm_name == "bfs"){
            bfs.m_enabled = true;
            bfs.m_source_vertex = stoull( props.get_property("bfs.source-vertex") );
        } else if (algorithm_name == "cdlp"){
            cdlp.m_enabled = true;
            cdlp.m_max_iterations = stoull( props.get_property("cdlp.max-iterations"));
        } else if (algorithm_name == "lcc"){
            lcc.m_enabled = true;
        } else if (algorithm_name == "pr") { // pagerank
            pagerank.m_enabled = true;
            pagerank.m_damping_factor = stod( props.get_property("pr.damping-factor"));
            pagerank.m_num_iterations = stoull( props.get_property("pr.num-iterations"));
        } else if(algorithm_name == "sssp"){
            sssp.m_enabled = true;
            sssp.m_source_vertex = stoull( props.get_property("sssp.source-vertex"));
        } else if(algorithm_name == "wcc"){
            wcc.m_enabled = true;
        }
    }
}

ostream& operator<<(std::ostream& out, const GraphalyticsAlgorithms& props){
    out << "[GraphalyticsAlgorithms";
    if(props.bfs.m_enabled){
        out << " BFS source: " << props.bfs.m_source_vertex << ";";
    }
    if(props.cdlp.m_enabled){
        out << " CDLP max_iterations: " << props.cdlp.m_max_iterations << ";";
    }
    if(props.lcc.m_enabled){
        out << " LCC;";
    }
    if(props.pagerank.m_enabled){
        out << " PageRank df: " << props.pagerank.m_damping_factor << ", num_iterations: " << props.pagerank.m_num_iterations << ";";
    }
    if(props.sssp.m_enabled){
        out << " SSSP source: " << props.sssp.m_source_vertex << ";";
    }
    if(props.wcc.m_enabled){
        out << " WCC;";
    }
    out << "]";
    return out;
}

/*****************************************************************************
 *                                                                           *
 *  GraphalyticsSequential                                                   *
 *                                                                           *
 ****************************************************************************/
GraphalyticsSequential::GraphalyticsSequential(std::shared_ptr<gfe::library::GraphalyticsInterface> interface, uint64_t num_repetitions, const GraphalyticsAlgorithms& properties) :
        m_interface(interface), m_num_repetitions(num_repetitions), m_properties(properties) { }

std::chrono::microseconds GraphalyticsSequential::execute(){
    auto interface = m_interface.get();
    constexpr uint64_t max_num_errors = 10; // if validation is enabled

    Timer t_global, t_local;
    t_global.start();

    for(uint64_t i = 0; i < m_num_repetitions; i++){

        if(m_properties.bfs.m_enabled){
            LOG("Execution " << (i+1) << "/" << m_num_repetitions << ": BFS from source vertex: " << m_properties.bfs.m_source_vertex);
            string path_tmp = get_temporary_path("bfs", i);
            const char* path_result = m_validate_output_enabled ? path_tmp.c_str() : nullptr;
            try {
                t_local.start();
                interface->bfs(m_properties.bfs.m_source_vertex, path_result);
                t_local.stop();
                LOG(">> BFS Execution time: " << t_local);
                m_exec_bfs.push_back(t_local.microseconds());

                if(m_validate_output_enabled){
                    string path_reference = get_validation_path("BFS");
                    if(common::filesystem::exists(path_reference)){
                        GraphalyticsValidate::bfs(path_tmp, path_reference, max_num_errors, get_validation_map());
                        LOG(">> Validation succeeded");
                        m_validate_results.emplace_back("bfs", ValidationResult::SUCCEEDED);
                        std::filesystem::remove(path_tmp);
                    } else if (i == 0) { // report it only the first time
                        LOG(">> Validation skipped, the reference file `" << path_reference << "' does not exist");
                        m_validate_results.emplace_back("bfs", ValidationResult::SKIPPED);
                    }
                }
            } catch (library::TimeoutError& e){
                LOG(">> BFS TIMEOUT");
                m_exec_bfs.push_back(-1);
                m_properties.bfs.m_enabled = false;
            } catch (utility::GraphalyticsValidateError& e){
            	LOG(">> Validation failed: " << e.what());
            	m_validate_results.emplace_back("bfs", ValidationResult::FAILED);
            	m_properties.bfs.m_enabled = false;
            }
        }

        if(m_properties.cdlp.m_enabled){
            LOG("Execution " << (i+1) << "/" << m_num_repetitions << ": CDLP, max_iterations: " << m_properties.cdlp.m_max_iterations);
            string path_tmp = get_temporary_path("cdlp", i);
            const char* path_result = m_validate_output_enabled ? path_tmp.c_str() : nullptr;
            try {
                t_local.start();
                interface->cdlp(m_properties.cdlp.m_max_iterations, path_result);
                t_local.stop();
                LOG(">> CDLP Execution time: " << t_local);
                m_exec_cdlp.push_back(t_local.microseconds());

                if(m_validate_output_enabled){
                    string path_reference = get_validation_path("CDLP");
                    if(common::filesystem::exists(path_reference)){
                        GraphalyticsValidate::cdlp(path_tmp, path_reference, max_num_errors, get_validation_map());
                        LOG(">> Validation succeeded");
                        m_validate_results.emplace_back("cdlp", ValidationResult::SUCCEEDED);
                        std::filesystem::remove(path_tmp);
                    } else if (i == 0) { // report it only the first time
                        LOG(">> Validation skipped, the reference file `" << path_reference << "' does not exist");
                        m_validate_results.emplace_back("cdlp", ValidationResult::SKIPPED);
                    }
                }
            } catch(library::TimeoutError& e){
                LOG(">> CDLP TIMEOUT");
                m_exec_cdlp.push_back(-1);
                m_properties.cdlp.m_enabled = false;
            } catch(utility::GraphalyticsValidateError& e){
                LOG(">> Validation failed: " << e.what());
                m_validate_results.emplace_back("cdlp", ValidationResult::FAILED);
                m_properties.cdlp.m_enabled = false;
            }
        }

        if(m_properties.lcc.m_enabled){
            LOG("Execution " << (i+1) << "/" << m_num_repetitions << ": LCC");
            string path_tmp = get_temporary_path("lcc", i);
            const char* path_result = m_validate_output_enabled ? path_tmp.c_str() : nullptr;
            try {
                t_local.start();
                interface->lcc(path_result);
                t_local.stop();
                LOG(">> LCC Execution time: " << t_local);
                m_exec_lcc.push_back(t_local.microseconds());

                if(m_validate_output_enabled){
                    string path_reference = get_validation_path("LCC");
                    if(common::filesystem::exists(path_reference)){
                        GraphalyticsValidate::lcc(path_tmp, path_reference, max_num_errors, get_validation_map());
                        LOG(">> Validation succeeded");
                        m_validate_results.emplace_back("lcc", ValidationResult::SUCCEEDED);
                        std::filesystem::remove(path_tmp);
                    } else if (i == 0) { // report it only the first time
                        LOG(">> Validation skipped, the reference file `" << path_reference << "' does not exist");
                        m_validate_results.emplace_back("lcc", ValidationResult::SKIPPED);
                    }
                }
            } catch(library::TimeoutError& e){
                LOG(">> LCC TIMEOUT");
                m_exec_lcc.push_back(-1);
                m_properties.lcc.m_enabled = false;
            } catch(utility::GraphalyticsValidateError& e){
                LOG(">> Validation failed: " << e.what());
                m_validate_results.emplace_back("lcc", ValidationResult::FAILED);
                m_properties.lcc.m_enabled = false;
            }
        }

        if(m_properties.pagerank.m_enabled){
            LOG("Execution " << (i+1) << "/" << m_num_repetitions << ": PageRank, damping factor: " << m_properties.pagerank.m_damping_factor << ", num_iterations: " << m_properties.pagerank.m_num_iterations);
            string path_tmp = get_temporary_path("pagerank", i);
            const char* path_result = m_validate_output_enabled ? path_tmp.c_str() : nullptr;
            try {
                t_local.start();
                interface->pagerank(m_properties.pagerank.m_num_iterations, m_properties.pagerank.m_damping_factor, path_result);
                t_local.stop();
                LOG(">> PageRank Execution time: " << t_local);
                m_exec_pagerank.push_back(t_local.microseconds());

                if(m_validate_output_enabled){
                    string path_reference = get_validation_path("PR");
                    if(common::filesystem::exists(path_reference)){
                        GraphalyticsValidate::pagerank(path_tmp, path_reference, max_num_errors, get_validation_map());
                        LOG(">> Validation succeeded");
                        m_validate_results.emplace_back("pagerank", ValidationResult::SUCCEEDED);
                        std::filesystem::remove(path_tmp);
                    } else if (i == 0) { // report it only the first time
                        LOG(">> Validation skipped, the reference file `" << path_reference << "' does not exist");
                        m_validate_results.emplace_back("pagerank", ValidationResult::SKIPPED);
                    }
                }
            } catch(library::TimeoutError& e){
                LOG(">> PageRank TIMEOUT");
                m_exec_pagerank.push_back(-1);
                m_properties.pagerank.m_enabled = false;
            } catch(utility::GraphalyticsValidateError& e){
                LOG(">> Validation failed: " << e.what());
                m_validate_results.emplace_back("pagerank", ValidationResult::FAILED);
                m_properties.pagerank.m_enabled = false;
            }
        }

        if(m_properties.sssp.m_enabled){
            LOG("Execution " << (i+1) << "/" << m_num_repetitions << ": SSSP, source: " << m_properties.sssp.m_source_vertex);
            string path_tmp = get_temporary_path("sssp", i);
            const char* path_result = m_validate_output_enabled ? path_tmp.c_str() : nullptr;
            try {
                t_local.start();
                interface->sssp(m_properties.sssp.m_source_vertex, path_result);
                t_local.stop();
                LOG(">> SSSP Execution time: " << t_local);
                m_exec_sssp.push_back(t_local.microseconds());

                if(m_validate_output_enabled){
                    string path_reference = get_validation_path("SSSP");
                    if(common::filesystem::exists(path_reference)){
                        GraphalyticsValidate::sssp(path_tmp, path_reference, max_num_errors, get_validation_map());
                        LOG(">> Validation succeeded");
                        m_validate_results.emplace_back("sssp", ValidationResult::SUCCEEDED);
                        std::filesystem::remove(path_tmp);
                    } else if (i == 0) { // report it only the first time
                        LOG(">> Validation skipped, the reference file `" << path_reference << "' does not exist");
                        m_validate_results.emplace_back("sssp", ValidationResult::SKIPPED);
                    }
                }
            } catch(library::TimeoutError& e){
                LOG(">> SSSP TIMEOUT");
                m_exec_sssp.push_back(-1);
                m_properties.sssp.m_enabled = false;
            } catch(utility::GraphalyticsValidateError& e){
                LOG(">> Validation failed: " << e.what());
                m_validate_results.emplace_back("sssp", ValidationResult::FAILED);
                m_properties.sssp.m_enabled = false;
            }
        }
        if(m_properties.wcc.m_enabled){
            LOG("Execution " << (i+1) << "/" << m_num_repetitions << ": WCC");
            string path_tmp = get_temporary_path("wcc", i);
            const char* path_result = m_validate_output_enabled ? path_tmp.c_str() : nullptr;
            try {
                t_local.start();
                interface->wcc(path_result);
                t_local.stop();
                LOG(">> WCC Execution time: " << t_local);
                m_exec_wcc.push_back(t_local.microseconds());

                if(m_validate_output_enabled){
                    string path_reference = get_validation_path("WCC");
                    if(common::filesystem::exists(path_reference)){
                        GraphalyticsValidate::wcc(path_tmp, path_reference, max_num_errors, get_validation_map());
                        LOG(">> Validation succeeded");
                        m_validate_results.emplace_back("wcc", ValidationResult::SUCCEEDED);
                        std::filesystem::remove(path_tmp);
                    } else if (i == 0) { // report it only the first time
                        LOG(">> Validation skipped, the reference file `" << path_reference << "' does not exist");
                        m_validate_results.emplace_back("wcc", ValidationResult::SKIPPED);
                    }
                }
            } catch(library::TimeoutError& e){
                LOG(">> WCC TIMEOUT");
                m_exec_wcc.push_back(-1);
                m_properties.wcc.m_enabled = false;
            } catch(utility::GraphalyticsValidateError& e){
                LOG(">> Validation failed: " << e.what());
                m_validate_results.emplace_back("wcc", ValidationResult::FAILED);
                m_properties.wcc.m_enabled = false;
            }
        }
    }

    t_global.stop();

    return t_global.duration<chrono::microseconds>();
}

void GraphalyticsSequential::report(bool save_in_db){
    if(!m_exec_bfs.empty()){
        ExecStatistics stats { m_exec_bfs };
        cout << ">> BFS " << stats << "\n";
        if(save_in_db) stats.save("bfs");
    }
    if(!m_exec_cdlp.empty()){
        ExecStatistics stats { m_exec_cdlp };
        cout << ">> CDLP " << stats << "\n";
        if(save_in_db) stats.save("cdlp");
    }
    if(!m_exec_lcc.empty()){
        ExecStatistics stats { m_exec_lcc };
        cout << ">> LCC " << stats << "\n";
        if(save_in_db) stats.save("lcc");
    }
    if(!m_exec_pagerank.empty()){
        ExecStatistics stats { m_exec_pagerank };
        cout << ">> PageRank " << stats << "\n";
        if(save_in_db) stats.save("pagerank");
    }
    if(!m_exec_sssp.empty()){
        ExecStatistics stats { m_exec_sssp };
        cout << ">> SSSP " << stats << "\n";
        if(save_in_db) stats.save("sssp");
    }
    if(!m_exec_wcc.empty()){
        ExecStatistics stats { m_exec_wcc };
        cout << ">> WCC " << stats << "\n";
        if(save_in_db) stats.save("wcc");
    }

    if(!m_validate_results.empty()){
        uint64_t num_validation_errors = 0;

    	for(const auto& pair : m_validate_results){
    	    num_validation_errors += (pair.second == ValidationResult::FAILED);

    	    if(save_in_db){
                auto store = configuration().db()->add("graphalytics_validation");
                store.add("algorithm", pair.first);
                switch(pair.second){
                case ValidationResult::SUCCEEDED:
                    store.add("result", "success");
                    break;
                case ValidationResult::FAILED:
                    store.add("result", "failure");
                    break;
                case ValidationResult::SKIPPED:
                    store.add("result", "skipped");
                    break;
                }
    	    }
    	}

    	if(num_validation_errors == 0){
    	    cout << ">> All executions succeeded the validation (" << m_validate_results.size() << " validation runs)" << endl;
    	} else {
    	    cout << ">> Executions that failed the validation: " << num_validation_errors << " out of " << m_validate_results.size() << " validation runs" << endl;
    	}
    }
}

void GraphalyticsSequential::set_validate_output(const std::string& path_property_file){
    string path = path_property_file;
    if(!common::filesystem::file_exists(path)){
        path += ".properties"; // try again by adding the suffix .properties
        if(!common::filesystem::file_exists(path)){
            ERROR("The file `" + path_property_file + "' does not exist");
        }
    }

    { // generate the temporary folder
        stringstream ss;
        ss << P_tmpdir << "/" << "graphalytics.validate_output.XXXXXX";
        string dynpath = ss.str();
        constexpr uint64_t buffer_sz = 1024;
        char buffer[buffer_sz];
        if(buffer_sz < dynpath.size() +1 /* \0 */ ) ERROR("Cannot generate the temporary path, the path name is too long");
        strcpy(buffer, dynpath.c_str());
        char* result = mkdtemp(buffer);
        if(result == nullptr){ ERROR("Cannot create the temporary path, mkdtemp error"); }
        m_validate_output_temp_dir = result;
        LOG("Temporary directory to store the output files: " << m_validate_output_temp_dir);
    }

    m_validate_path_expected = path;
    m_validate_output_enabled = true;
}

void GraphalyticsSequential::set_validate_remap_vertices(const std::string& path_property_file){
    string path_results = path_property_file;
    if(!common::filesystem::file_exists(path_results)){
        path_results += ".properties"; // try again by adding the suffix .properties
        if(!common::filesystem::file_exists(path_results)){
            ERROR("The file `" + path_property_file + "' does not exist");
        }
    }
    Timer timer; timer.start();
    LOG("Validation: mapping the vertices from `" << m_validate_path_expected << "' into `" << path_results << "' ... ");

    reader::GraphalyticsReader reader_expected { m_validate_path_expected };
    reader::GraphalyticsReader reader_results { path_results };

    uint64_t expected_num_vertices = stoull( reader_expected.get_property("meta.vertices") );
    uint64_t results_num_vertices = stoull( reader_results.get_property("meta.vertices") );

    if(expected_num_vertices != results_num_vertices) {
        ERROR("The number of vertices in the two files `" << m_validate_path_expected << "' and `" << path_results << "' do not match: "
                << expected_num_vertices << " != " << results_num_vertices)
    }

    // the vertices in the results file must have been remapped following the same sorted order of the expected file
    bool vtx_stable_map = reader_results.get_property("meta.stable-map") == "true";
    if(!vtx_stable_map){
        ERROR("[Validation] We cannot compare `" << m_validate_path_expected << "' and `" << path_results << "' because the vertices have "
               "not been mapped following the same sorted order of the expected input graph. The property `meta.stable-map' is not set or it is false "
               "in the file: `" << path_results << "'.")
    }

    m_validate_path_expected.reserve(expected_num_vertices);

    uint64_t vertex_expected, vertex_results;
    while(reader_expected.read_vertex(vertex_expected) && reader_results.read_vertex(vertex_results)){
        m_validation_map[vertex_expected] = vertex_results;
    }

    if(reader_expected.read_vertex(vertex_expected)){
        ERROR("The file `" << m_validate_path_expected << "' contains more vertices than `" << path_results << "'");
    }
    if(reader_results.read_vertex(vertex_results)){
        ERROR("The file `" << m_validate_path_expected << "' contains less vertices than `" << path_results << "'");
    }

    timer.stop();
    LOG("Validation: mapping completed in " << timer);
}

string GraphalyticsSequential::get_temporary_path(const string& algorithm_name, uint64_t execution_no) const{
    if(!m_validate_output_enabled) return "";

    stringstream ss;
    ss << m_validate_output_temp_dir << "/";
    ss << algorithm_name << "_" << execution_no << ".txt";
    return ss.str();
}

string GraphalyticsSequential::get_validation_path(const string& algorithm_suffix) const {
    if(!m_validate_output_enabled) return "";
    string path = m_validate_path_expected;
    auto index_of = path.find_last_of('.');
    if(index_of != string::npos && path.substr(index_of) == ".properties"){
        path = path.substr(0, index_of);
    }
    return path + "-" + algorithm_suffix;
}

const std::unordered_map<uint64_t, uint64_t>* GraphalyticsSequential::get_validation_map() const {
    if(m_validation_map.empty()){ // no mapping required
        return nullptr;
    } else {
        return &m_validation_map;
    }
}

} // namespace experiment
