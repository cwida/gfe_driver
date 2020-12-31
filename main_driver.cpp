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
#include <iostream>

#include "common/database.hpp"
#include "common/error.hpp"
#include "common/system.hpp"
#include "common/timer.hpp"
#include "experiment/aging2_experiment.hpp"
#include "experiment/insert_only.hpp"
#include "experiment/graphalytics.hpp"
#include "experiment/validate.hpp"
#include "graph/edge_stream.hpp"
#include "library/interface.hpp"
#include "third-party/cxxopts/cxxopts.hpp"
#include "utility/memory_usage.hpp"

#include "configuration.hpp"
#if defined(HAVE_OPENMP)
#include "omp.h"
#endif

extern "C" { void invoke_me(); }

using namespace common;
using namespace gfe;
using namespace gfe::experiment;
using namespace std;

static void run_standalone(int argc, char* argv[]){
    configuration().initialise(argc, argv);
    if(configuration().get_aging_memfp() && !configuration().get_aging_memfp_physical()){
        utility::MemoryUsage::initialise(argc, argv); // init the memory profiler
    }
    LOG("[driver] Initialising ...");

    if(configuration().has_database()){
        LOG( "[driver] Save the current configuration properties in " << configuration().get_database_path() )
        configuration().save_parameters();
    }

    const std::string& path_graph = configuration().get_path_graph(); // graph to use?
    if(path_graph.empty()) ERROR("Path to the graph to load not set (use the parameter --graph)");

    // OpenMP settings
#if defined(HAVE_OPENMP)
    if(configuration().num_threads_omp() > 0){
        LOG("[driver] OpenMP, max number of threads: " << configuration().num_threads_omp());
        omp_set_num_threads(configuration().num_threads_omp());
    }
#endif

    // implementation to evaluate
    LOG("[driver] Library name: " << configuration().get_library_name() );
    shared_ptr<library::Interface> impl { configuration().generate_graph_library() };
    impl->set_timeout(configuration().get_timeout_graphalytics());


    auto impl_ga = dynamic_pointer_cast<library::GraphalyticsInterface>(impl);
    if(impl_ga.get() == nullptr && configuration().num_repetitions() > 0){ // Shall we execute the Graphalytics suite?
        ERROR("The library does not support the Graphalytics suite of algorithms");
    }

    LOG("[driver] The library is set for a directed graph: " << (configuration().is_graph_directed() ? "yes" : "no"));

    uint64_t random_vertex = numeric_limits<uint64_t>::max();
    int64_t num_validation_errors = -1; // -1 => no validation performed
    if(configuration().is_load()){
        auto impl_load = dynamic_pointer_cast<library::LoaderInterface>(impl);
        if(impl_load.get() == nullptr){ ERROR("The library `" << configuration().get_library_name() << "' does not support loading"); }
        auto impl_rndvtx = dynamic_pointer_cast<library::RandomVertexInterface>(impl);
        if(impl_rndvtx.get() == nullptr){ ERROR("The library `" << configuration().get_library_name() << "' does not allow to fetch a random vertex"); }

        LOG("[driver] Loading the graph: " << path_graph);
        common::Timer timer; timer.start();
        impl_load->load(path_graph);
        timer.stop();
        LOG("[driver] Load performed in " << timer);

        if(configuration().validate_inserts() && impl_load->can_be_validated()){
            auto stream = make_shared<graph::WeightedEdgeStream> ( configuration().get_path_graph() );
            num_validation_errors = validate_updates(impl_load, stream);
        }

        random_vertex = impl_rndvtx->get_random_vertex_id();
    } else {
        auto impl_upd = dynamic_pointer_cast<library::UpdateInterface>(impl);
        if(impl_upd.get() == nullptr){ ERROR("The library `" << configuration().get_library_name() << "' does not support updates"); }

        if(configuration().get_update_log().empty()){
            LOG("[driver] Using the graph " << path_graph);
            auto stream = make_shared<graph::WeightedEdgeStream> ( configuration().get_path_graph() );
            stream->permute();
            if(stream->num_edges() > 0) random_vertex = stream->get(0).m_source;

            LOG("[driver] Number of concurrent threads: " << configuration().num_threads(THREADS_WRITE) );

            if(configuration().measure_latency()) ERROR("[driver] InsertOnly, support for latency measurements removed");

            InsertOnly experiment { impl_upd, stream, configuration().num_threads(THREADS_WRITE) };
            experiment.set_build_frequency(chrono::milliseconds{ configuration().get_build_frequency() });
            experiment.set_scheduler_granularity(1ull < 20);
            experiment.execute();
            if(configuration().has_database()) experiment.save();

            if(configuration().validate_inserts() && impl_upd->can_be_validated()){
                num_validation_errors = validate_updates(impl_upd, stream);
            }

        } else {
            LOG("[driver] Number of concurrent threads: " << configuration().num_threads(THREADS_WRITE) );
            LOG("[driver] Aging2, path to the log of updates: " << configuration().get_update_log());
            Aging2Experiment experiment;
            experiment.set_library(impl_upd);
            experiment.set_log(configuration().get_update_log());
            experiment.set_parallelism_degree(configuration().num_threads(THREADS_WRITE));
            experiment.set_release_memory(configuration().get_aging_release_memory());
            experiment.set_report_progress(true);
            experiment.set_report_memory_footprint(configuration().get_aging_memfp_report());
            experiment.set_build_frequency(chrono::milliseconds{ configuration().get_build_frequency() });
            experiment.set_max_weight(configuration().max_weight());
            experiment.set_measure_latency(configuration().measure_latency());
            experiment.set_num_reports_per_ops(configuration().get_num_recordings_per_ops());
            experiment.set_timeout(chrono::seconds { configuration().get_timeout_aging2() });
            experiment.set_measure_memfp( configuration().measure_memfp() );
            experiment.set_memfp_physical( configuration().get_aging_memfp_physical() );
            experiment.set_memfp_threshold(configuration().get_aging_memfp_threshold());
            experiment.set_cooloff(chrono::seconds { configuration().get_aging_cooloff_seconds() });

            auto result = experiment.execute();
            if(configuration().has_database()) result.save(configuration().db());
            random_vertex = result.get_random_vertex_id();

            if(configuration().validate_inserts() && impl_upd->can_be_validated()){
                LOG("[driver] Validation of updates requested, loading the original graph from: " << path_graph);
                auto stream = make_shared<graph::WeightedEdgeStream> ( configuration().get_path_graph() );
                num_validation_errors = validate_updates(impl_upd, stream);
            }
        }
    }

    if(configuration().has_database()){
        vector<pair<string, string>> params;
        params.push_back(make_pair("num_validation_errors", to_string(num_validation_errors)));
        configuration().db()->store_parameters(params);
    }

    if(configuration().num_repetitions() > 0){

#if defined(HAVE_OPENMP)
        if(configuration().num_threads(ThreadsType::THREADS_READ) != 0 ){
            LOG("[driver] OpenMP, number of threads for the Graphalytics suite: " << configuration().num_threads(ThreadsType::THREADS_READ));
            omp_set_num_threads(configuration().num_threads(ThreadsType::THREADS_READ));
        }
#endif

        // run the graphalytics suite
        GraphalyticsAlgorithms properties { path_graph };

        if(properties.bfs.m_enabled == true && properties.sssp.m_enabled == false){
            LOG("[driver] Enabling SSSP with random weights, source vertex: " << random_vertex);
            properties.sssp.m_enabled = true;
            properties.sssp.m_source_vertex = random_vertex;
        }

        configuration().blacklist(properties);

        GraphalyticsSequential exp_seq { impl_ga, configuration().num_repetitions(), properties };

        if(configuration().validate_output()){
            LOG("[driver] Enabling validation mode");
            exp_seq.set_validate_output( configuration().get_validation_graph() );
            if(configuration().get_validation_graph() != path_graph){
                exp_seq.set_validate_remap_vertices( path_graph );
            }
        }

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

    return rc;
}


