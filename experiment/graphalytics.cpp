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
#include <iostream>
#include <string>
#include <sstream>

#include "common/filesystem.hpp"
#include "common/timer.hpp"
#include "library/interface.hpp"
#include "reader/graphalytics_reader.hpp"
#include "configuration.hpp"
#include "statistics.hpp"

using namespace common;
using namespace std;

namespace experiment {

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
GraphalyticsAlgorithms::GraphalyticsAlgorithms(const string& path){
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
GraphalyticsSequential::GraphalyticsSequential(std::shared_ptr<library::GraphalyticsInterface> interface, uint64_t num_repetitions, const GraphalyticsAlgorithms& properties) :
        m_interface(interface), m_num_repetitions(num_repetitions), m_properties(properties) { }

std::chrono::microseconds GraphalyticsSequential::execute(){
    auto interface = m_interface.get();

    Timer t_global, t_local;
    t_global.start();

    for(uint64_t i = 0; i < m_num_repetitions; i++){
        if(m_properties.bfs.m_enabled){
            LOG("Execution " << (i+1) << "/" << m_num_repetitions << ": BFS from source vertex: " << m_properties.bfs.m_source_vertex);
            t_local.start();
            interface->bfs(m_properties.bfs.m_source_vertex);
            t_local.stop();
            LOG(">> BFS Execution time: " << t_local);
            m_exec_bfs.push_back(t_local.microseconds());
        }
        if(m_properties.cdlp.m_enabled){
            LOG("Execution " << (i+1) << "/" << m_num_repetitions << ": CDLP, max_iterations: " << m_properties.cdlp.m_max_iterations);
            t_local.start();
            interface->cdlp(m_properties.cdlp.m_max_iterations);
            t_local.stop();
            LOG(">> CDLP Execution time: " << t_local);
            m_exec_cdlp.push_back(t_local.microseconds());
        }
        if(m_properties.lcc.m_enabled){
            LOG("Execution " << (i+1) << "/" << m_num_repetitions << ": LCC");
            t_local.start();
            interface->lcc();
            t_local.stop();
            LOG(">> LCC Execution time: " << t_local);
            m_exec_lcc.push_back(t_local.microseconds());
        }
        if(m_properties.pagerank.m_enabled){
            LOG("Execution " << (i+1) << "/" << m_num_repetitions << ": PageRank, damping factor: " << m_properties.pagerank.m_damping_factor << ", num_iterations: " << m_properties.pagerank.m_num_iterations);
            t_local.start();
            interface->pagerank(m_properties.pagerank.m_num_iterations, m_properties.pagerank.m_damping_factor);
            t_local.stop();
            LOG(">> PageRank Execution time: " << t_local);
            m_exec_pagerank.push_back(t_local.microseconds());
        }
        if(m_properties.sssp.m_enabled){
            LOG("Execution " << (i+1) << "/" << m_num_repetitions << ": SSSP, source: " << m_properties.sssp.m_source_vertex);
            t_local.start();
            interface->sssp(m_properties.sssp.m_source_vertex);
            t_local.stop();
            LOG(">> SSSP Execution time: " << t_local);
            m_exec_sssp.push_back(t_local.microseconds());
        }
        if(m_properties.wcc.m_enabled){
            LOG("Execution " << (i+1) << "/" << m_num_repetitions << ": WCC");
            t_local.start();
            interface->wcc();
            t_local.stop();
            LOG(">> WCC Execution time: " << t_local);
            m_exec_wcc.push_back(t_local.microseconds());
        }
    }

    t_global.stop();

    return t_global.duration<chrono::microseconds>();
}

void GraphalyticsSequential::report(bool save_in_db){
    if(m_properties.bfs.m_enabled){
        ExecStatistics stats { m_exec_bfs };
        cout << ">> BFS " << stats << "\n";
        if(save_in_db) stats.save("bfs");
    }
    if(m_properties.cdlp.m_enabled){
        ExecStatistics stats { m_exec_cdlp };
        cout << ">> CDLP " << stats << "\n";
        if(save_in_db) stats.save("cdlp");
    }
    if(m_properties.lcc.m_enabled){
        ExecStatistics stats { m_exec_lcc };
        cout << ">> LCC " << stats << "\n";
        if(save_in_db) stats.save("lcc");
    }
    if(m_properties.pagerank.m_enabled){
        ExecStatistics stats { m_exec_pagerank };
        cout << ">> PageRank " << stats << "\n";
        if(save_in_db) stats.save("pagerank");
    }
    if(m_properties.sssp.m_enabled){
        ExecStatistics stats { m_exec_sssp };
        cout << ">> SSSP " << stats << "\n";
        if(save_in_db) stats.save("sssp");
    }
    if(m_properties.wcc.m_enabled){
        ExecStatistics stats { m_exec_wcc };
        cout << ">> WCC " << stats << "\n";
        if(save_in_db) stats.save("wcc");
    }
}

} // namespace experiment
