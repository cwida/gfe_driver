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

//#include <algorithm>
//#include <cctype>
//#include <cstdlib>
#include <iostream>
//
#include "common/error.hpp"
//#include "common/timer.hpp"
#include "experiment/aging.hpp"
#include "experiment/insert_only.hpp"
#include "experiment/graphalytics.hpp"
#include "graph/edge_stream.hpp"
#include "library/interface.hpp"
#include "third-party/cxxopts/cxxopts.hpp"

#include "configuration.hpp"

using namespace common;
using namespace experiment;
using namespace std;

static void run_standalone(int argc, char* argv[]){
    LOG("[driver] Init configuration ... " );
    StandaloneConfiguration::initialise(argc, argv);

    if(configuration().has_database()){
        LOG( "[driver] Save the current configuration properties in " << configuration().get_database_path() )
        configuration().save_parameters();
    }

    const std::string& path_graph = cfgdriver().get_path_graph(); // graph to use?
    if(path_graph.empty()) ERROR("Path to the graph to load not set (use the parameter --graph)");

    // implementation to the evaluate
    LOG("[driver] Library name: " << cfgstandalone().get_library_name() );
    shared_ptr<library::Interface> impl { cfgstandalone().generate_graph_library() };
    impl->set_timeout(cfgdriver().get_timeout_per_operation());

    auto impl_upd = dynamic_pointer_cast<library::UpdateInterface>(impl);
    if(impl_upd.get() == nullptr){ ERROR("The library does not support updates"); }

    auto impl_ga = dynamic_pointer_cast<library::GraphalyticsInterface>(impl);
    if(impl_ga.get() == nullptr && cfgdriver().num_repetitions() > 0){ // Shall we execute the Graphalytics suite?
        ERROR("The library does not support the Graphalytics suite of algorithms");
    }

    LOG("[driver] The library is set for a directed graph: " << (cfgstandalone().is_graph_directed() ? "yes" : "no"));

    LOG("[driver] Loading the graph from " << path_graph);
    auto stream = make_shared<graph::WeightedEdgeStream> ( cfgdriver().get_path_graph() );
    stream->permute();
    uint64_t random_vertex = stream->num_edges() > 0 ? stream->get(0).m_source : 0;

    LOG("[driver] Number of concurrent threads: " << cfgdriver().num_threads(THREADS_WRITE) );
    if(cfgdriver().coefficient_aging() == 0.0){ // insert the elements in the graph one by one
        InsertOnly experiment { impl_upd, move(stream), cfgdriver().num_threads(THREADS_WRITE) };
        experiment.execute();
        if(configuration().has_database()) experiment.save();
    } else {
        LOG("[driver] Number of updates to perform: " << stream->num_edges() * cfgdriver().coefficient_aging());
        Aging experiment(impl_upd, move(stream), cfgdriver().coefficient_aging(), cfgdriver().num_threads(THREADS_WRITE));
        experiment.set_expansion_factor_vertices(cfgdriver().get_ef_vertices());
        experiment.set_expansion_factor_edges(cfgdriver().get_ef_edges());
        experiment.set_report_progress(true);
        experiment.set_build_frequency(chrono::milliseconds{ cfgdriver().get_build_frequency() });
        experiment.execute();
        if(configuration().has_database()) experiment.save();
    }

    if(cfgdriver().num_repetitions() > 0){
        // run the graphalytics suite
        GraphalyticsAlgorithms properties { path_graph };

        if(properties.bfs.m_enabled == true && properties.sssp.m_enabled == false){
            LOG("[driver] Enabling SSSP with random weights, source vertex: " << random_vertex);
            properties.sssp.m_enabled = true;
            properties.sssp.m_source_vertex = random_vertex;
        }

        GraphalyticsSequential exp_seq { impl_ga, cfgdriver().num_repetitions(), properties };
        exp_seq.execute();
        exp_seq.report(configuration().has_database());
    }

    LOG( "[driver] Done" );
}

int main(int argc, char* argv[]){
    int rc = 0;
    try {
        run_standalone(argc, argv);
    } catch(common::Error& e){
        cerr << e << endl;
        cerr << "Client terminating due to exception..." << endl;
        rc = 1;
    } catch(cxxopts::option_not_exists_exception& e){
        cerr << "ERROR: " << e.what() << "\n";
        rc = 1;
    } catch(cxxopts::option_requires_argument_exception& e){
        cerr << "ERROR: Invalid command line option, " << e.what() << "\n";
        rc = 1;
    }

    cfgfree();
    return rc;
}


