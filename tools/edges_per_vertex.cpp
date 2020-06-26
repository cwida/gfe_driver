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

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <memory>
#include <omp.h>
#include <shared_mutex>
#include <sstream>
#include <vector>

// libcommon
#include "common/cpu_topology.hpp"
#include "common/error.hpp"
#include "common/filesystem.hpp"
#include "common/permutation.hpp"
#include "common/quantity.hpp"
#include "common/system.hpp"
#include "common/timer.hpp"

// gfe
#include "experiment/insert_only.hpp"
#include "graph/edge.hpp"
#include "graph/edge_stream.hpp"
#include "graph/vertex_list.hpp"
#include "library/interface.hpp"
#include "reader/graphalytics_reader.hpp"
#include "configuration.hpp"

// teseo
#if defined(HAVE_TESEO)
#include "library/teseo/teseo_driver.hpp"
#include "teseo/context/global_context.hpp"
#include "teseo/memstore/memstore.hpp"
#include "teseo.hpp"
#endif

// graphone
#if defined(HAVE_GRAPHONE)
#include "library/graphone/graphone.hpp"
#include "library/graphone/internal.hpp"
#endif

// llama
#if defined(HAVE_LLAMA)
#include "library/llama/llama_ref.hpp"
#include "library/llama/llama_internal.hpp"
#endif

#if defined(HAVE_STINGER)
#include "library/stinger/stinger.hpp"
#include "stinger_core/stinger.h"
#endif

using namespace gfe;
using namespace std;

// globals
static string g_destination;
static string g_library;
static string g_path_graph;
static shared_ptr<gfe::library::Interface> g_interface { nullptr };

// function prototypes
static void load();
static void parse_args(int argc, char* argv[]);
static void run();
[[maybe_unused]] static void run_teseo();
[[maybe_unused]] static void run_graphone();
[[maybe_unused]] static void run_llama();
[[maybe_unused]] static void run_stinger();
static string string_usage(char* program_name);
int main(int argc, char* argv[]){
    parse_args(argc, argv);
    load();
    cout << "Library: " << g_library << ", destination: " << g_destination << " ... " << endl;
    run();
    cout << "\nDone" << endl;
    return 0;
}

static void run(){
    common::Timer timer;
    timer.start();

    if(g_library == "teseo"){
#if defined(HAVE_TESEO)
        run_teseo();
#else
        assert(0 && "Support for teseo disabled");
#endif
    } else if(g_library == "graphone"){
#if defined(HAVE_GRAPHONE)
        run_graphone();
#else
        assert(0 && "Support for graphone disabled");
#endif
    } else if(g_library == "llama"){
#if defined(HAVE_LLAMA)
        run_llama();
#else
        assert(0 && "Support for llama disabled");
#endif
    } else if(g_library == "stinger"){
#if defined(HAVE_STINGER)
        run_stinger();
#else
        assert(0 && "Support for stinger disabled");
#endif
    } else {
        assert(0 && "Invalid library");
    }

    timer.stop();
    LOG("Execution completed in " << timer);
}

#if defined(HAVE_TESEO)
static void run_teseo(){
    fstream file (g_destination.c_str(), ios_base::out);
    if(!file.good()){ ERROR("Invalid file name: " << file); }

    using namespace teseo;
    Teseo* teseo = reinterpret_cast<Teseo*>(dynamic_cast<library::TeseoDriver*>(g_interface.get())->handle_impl());
    teseo->register_thread();
    auto tx = teseo->start_transaction(/* read only ? */ true);
    for(uint64_t i = 0, N = tx.num_vertices(); i < N; i++){
        file << tx.degree(i, /* logical ? */ true) << "\n";
    }
    tx.commit();
    teseo->unregister_thread();

    file.close();
}
#endif


#if defined(HAVE_LLAMA)
static void _bm_run_llama(){
    fstream file (g_destination.c_str(), ios_base::out);
    if(!file.good()){ ERROR("Invalid file name: " << file); }

    auto instance = dynamic_cast<library::LLAMAClass*>(g_interface.get());
    shared_lock<library::LLAMAClass::shared_mutex_t> slock(instance->m_lock_checkpoint);
    auto graph = instance->get_snapshot();

    for(uint64_t i = 0, N = instance->num_vertices(); i < N; i++){
        file << graph.out_degree(i) + graph.in_degree(i) << "\n";
    }

    file.close();
}

static void run_llama(){
    _bm_run_llama();
}
#endif

#if defined(HAVE_GRAPHONE)
static void run_graphone(){
    fstream file (g_destination.c_str(), ios_base::out);
    if(!file.good()){ ERROR("Invalid file name: " << file); }

    auto* view = create_static_view(get_graphone_graph(), SIMPLE_MASK | PRIVATE_MASK); // global
    uint64_t num_vertices = g_interface->num_vertices();

    for(uint64_t i = 0; i < num_vertices; i++){
        file << view->get_degree_out(i) << "\n";
    }

    delete_static_view(view);

    file.close();
}
#endif

#if defined(HAVE_STINGER)
static void run_stinger(){
    fstream file (g_destination.c_str(), ios_base::out);
    if(!file.good()){ ERROR("Invalid file name: " << file); }

    auto stinger = reinterpret_cast<struct stinger*>(dynamic_cast<library::StingerRef*>(g_interface.get())->handle());
    uint64_t num_vertices = g_interface->num_vertices();

    for(uint64_t i = 0; i < num_vertices; i++){
        file << stinger_outdegree_get(stinger, i) << "\n";
    }

    file.close();
}
#endif

static void load(){
    if(g_library == "teseo"){
#if defined(HAVE_TESEO)
        g_interface.reset( new gfe::library::TeseoDriver(/* directed ? */ false) );
#else
        cerr << "ERROR: gfe configured and built without linking the library teseo\n";
        exit(EXIT_FAILURE);
#endif
    } else if(g_library == "graphone"){
#if defined(HAVE_GRAPHONE)
        uint64_t ram = common::get_total_ram(); // in bytes
        // Ensure that the vertex array does not use more than 20% of the available RAM
        uint64_t num_vertices = std::min<uint64_t>(4ull<<30, ram/(5 * 8)); // 8 bytes per vertex, 1/5 of the total ram

        { // critical section, for the output
            std::scoped_lock<std::mutex> lock{_log_mutex};
            cout << "GraphOne: capacity of the vertex array implicitly set to: " << common::ComputerQuantity(num_vertices) << " vertices " << endl;
        }

        g_interface.reset( new gfe::library::GraphOne(/* directed ? */ false, /* vtx dict ? */ true, /* blind writes ? */ true, /* ignore build ? */ false, /* ref impl ? */ true, num_vertices) );
#else
        cerr << "ERROR: gfe configured and built without linking the library graphone\n";
        exit(EXIT_FAILURE);
#endif
    } else if(g_library == "llama"){
#if defined(HAVE_LLAMA)
        g_interface.reset( new library::LLAMARef(/* directed ? */ false) );
#else
        cerr << "ERROR: gfe configured and built without linking the library llama\n";
        exit(EXIT_FAILURE);
#endif
    } else if(g_library == "stinger"){
#if defined(HAVE_STINGER)
        g_interface.reset( new library::StingerRef(/* directed ? */ false) );
#else
        cerr << "ERROR: gfe configured and built without linking the library stinger\n";
        exit(EXIT_FAILURE);
#endif
    }

    LOG("Loading the graph from " << g_path_graph << " ...");
    auto edges = make_shared<gfe::graph::WeightedEdgeStream> ( g_path_graph );
    edges->permute();

    uint64_t num_threads = thread::hardware_concurrency();
    if(g_library == "stinger"){ // best number of threads in stones2 according to the scalability results
        num_threads = 1;
    } else if(g_library == "llama"){
        num_threads = 16;
    } else if(g_library == "graphone"){
        num_threads = 3;
    }

    LOG("Inserting " << edges->num_edges() << " edges into `" << g_library << "' ...");
    gfe::experiment::InsertOnly insert(dynamic_pointer_cast<gfe::library::UpdateInterface>(g_interface), edges, num_threads);
    if(g_library == "llama"){ insert.set_build_frequency( 10s ); }
    insert.set_scheduler_granularity(1ull < 20);
    insert.execute();
}

static void parse_args(int argc, char* argv[]) {
    static struct option long_options[] = {
        {"graph", required_argument, nullptr, 'G'},
        {"help", no_argument, nullptr, 'h'},
        {"library", required_argument, nullptr, 'l'},
        {0, 0, 0, 0} // keep at the end
    };

    int option { 0 };
    int option_index = 0;
    while( (option = getopt_long(argc, argv, "G:hl:", long_options, &option_index)) != -1 ){
        switch(option){
        case 'G': {
            string path_graph = optarg;
            if(!common::filesystem::file_exists(path_graph)){
                cerr << "ERROR: The file `" << path_graph << "' does not exist" << endl;
                exit(EXIT_FAILURE);
            }
            if(common::filesystem::extension(path_graph) != ".properties"){
                cerr << "ERROR: The file `" << path_graph << "' does not with the extension '.properties'. Only graphs from the Graphalytics set are supported." << endl;
                exit(EXIT_FAILURE);
            }
            g_path_graph = path_graph;
        } break;
        case 'h': {
            cout << "Set of micro benchmarks to evaluate point lookups and scans in Teseo\n";
            cout << string_usage(argv[0]) << endl;
            exit(EXIT_SUCCESS);
        } break;
        case 'l': {
            string library = optarg;
            if(library != "teseo" && library != "graphone" && library != "llama" && library != "stinger"){
                cerr << "ERROR: Invalid library: `" << library << "'. Only \"teseo\", \"graphone\", \"llama\" and \"stinger\" are supported." << endl;
            }
            g_library = library;
        } break;
        default:
            assert(0 && "Invalid option");
        }
    }

    if(optind < argc){
        g_destination = argv[optind];
    } else {
        cerr << "ERROR: output file not set\n";
        exit(EXIT_FAILURE);
    }

    if(g_path_graph.empty()){
        cerr << "ERROR: Input graph (-G) not specified\n";
        cerr << string_usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    if(g_library.empty()){
        g_library = "teseo";
    }
}

static string string_usage(char* program_name) {
    stringstream ss;
    ss << "Usage: " << program_name << " -G <graph> [-l <library>] <destination>\n";
    ss << "Where: \n";
    ss << "  -G <graph> is an .properties file of an undirected graph from the Graphalytics data set\n";
    ss << "  -l <library> is the library to execute. Only \"teseo\" (default), \"graphone\", \"llama\" and \"stinger\" are supported\n";
    ss << "  <destination> is the path where to store the output file\n";
    return ss.str();
}

