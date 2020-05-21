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
#include <sstream>
#include <vector>

// libcommon
#include "common/cpu_topology.hpp"
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

using namespace gfe;
using namespace std;

namespace {
struct Sample {
    string m_experiment; // the name of the experiment
    int m_num_threads; // parallelism degree
    uint64_t m_microsecs; // completion time, in microsecs

    Sample(const string& experiment, int num_threads, uint64_t microsecs) : m_experiment(experiment), m_num_threads(num_threads), m_microsecs(microsecs) { }
};
} // anon

// globals
static string g_library;
static int g_num_repetitions;
static vector<int> g_num_threads;
static string g_path_graph;
static uint64_t g_sum_degree;
static uint64_t g_sum_point_lookups;
static uint64_t g_sum_scan;
static vector<uint64_t> g_vertices_logical; // logical vertices, unsorted
static vector<uint64_t> g_vertices_sorted;
static vector<uint64_t> g_vertices_unsorted;
static shared_ptr<gfe::library::Interface> g_interface { nullptr };
static vector<Sample> g_samples;
static vector<string> g_experiments; // list of all experiments performed
static unordered_map</* experiment_name @ num_threads*/ string, /* median in microsecs */ uint64_t> g_medians;

// function prototypes
static void compute_medians(); // populate g_medians
static void load();
static void parse_args(int argc, char* argv[]);
static void run();
[[maybe_unused]] static void run_teseo();
[[maybe_unused]] static void run_graphone();
static void print_results();
static void save_results(const std::string& where);
static string string_usage(char* program_name);
template <typename Duration> static string to_string(Duration duration);
static void validate_sum_degree(uint64_t sum);
static void validate_sum_point_lookups(uint64_t sum);
static void validate_sum_scan(uint64_t sum);


int main(int argc, char* argv[]){
    parse_args(argc, argv);
    load();
    run();
    compute_medians();

    cout << "Library: " << g_library << "\n";

    print_results();

    string path_results = "/tmp/bm_";
    path_results += to_string(common::concurrency::get_process_id());
    path_results += ".json";
    save_results(path_results);
    cout << "Results saved into `" << path_results << "'\n";

    cout << "\nDone" << endl;
    return 0;
}

static void run(){
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
    } else {
        assert(0 && "Invalid library");
    }
}

#if defined(HAVE_TESEO)

namespace {

struct RegisterThread {
    RegisterThread& operator=(const RegisterThread&) = delete;
    teseo::Teseo* m_teseo;

public:
    RegisterThread(teseo::Teseo* teseo) : m_teseo(teseo){
        assert(teseo != nullptr);
        assert(omp_get_thread_num() == 0 && "Expected to be initialised in the master thread");
    }

    RegisterThread(const RegisterThread& rt) : m_teseo(rt.m_teseo){
        if(omp_get_thread_num() > 0){ m_teseo->register_thread(); }
    }

    ~RegisterThread(){
        if(omp_get_thread_num() > 0){ m_teseo->unregister_thread(); }
    }
};

} // anon

static void run_teseo(){
    using namespace teseo;
    Teseo* teseo = reinterpret_cast<Teseo*>(dynamic_cast<library::TeseoDriver*>(g_interface.get())->handle_impl());
    teseo->register_thread();
    auto tx_ro = teseo->start_transaction(/* read only ? */ true);
    auto iter_ro = tx_ro.iterator();
    const uint64_t num_vertices = g_vertices_sorted.size();
    common::Timer timer;
    RegisterThread rt(teseo);

    for(int r = 0; r < g_num_repetitions; r++){
        LOG("Repetition: " << (r +1) << "/" << g_num_repetitions);
        for(auto num_threads: g_num_threads){
            LOG("    num threads: " << num_threads);

            // vertices, sorted, vertex identifiers
            timer.start();
            uint64_t sum = 0;
            #pragma omp parallel for num_threads(num_threads) reduction(+:sum) firstprivate(rt, iter_ro)
            for(uint64_t i = 0; i < num_vertices; i++){
                sum += tx_ro.degree(g_vertices_sorted[i], false);
            }
            timer.stop();

            validate_sum_degree(sum);
            g_samples.emplace_back("degree_vtx_sorted", num_threads, timer.microseconds());

            // vertices, unsorted, vertex identifiers
            timer.start();
            sum = 0;
            #pragma omp parallel for num_threads(num_threads) reduction(+:sum) firstprivate(rt, iter_ro)
            for(uint64_t i = 0; i < num_vertices; i++){
                sum += tx_ro.degree(g_vertices_unsorted[i], false);
            }
            timer.stop();
            validate_sum_degree(sum);
            g_samples.emplace_back("degree_vtx_unsorted", num_threads, timer.microseconds());

            // vertices, logical identifiers, sorted
            timer.start();
            sum = 0;
            #pragma omp parallel for num_threads(num_threads) reduction(+:sum) firstprivate(rt, iter_ro)
            for(uint64_t i = 0; i < num_vertices; i++){
                sum += tx_ro.degree(i, true);
            }
            timer.stop();
            validate_sum_degree(sum);
            g_samples.emplace_back("degree_logical_sorted", num_threads, timer.microseconds());

            // vertices, logical identifiers, unsorted
            timer.start();
            sum = 0;
            #pragma omp parallel for num_threads(num_threads) reduction(+:sum) firstprivate(rt, iter_ro)
            for(uint64_t i = 0; i < num_vertices; i++){
                sum += tx_ro.degree(g_vertices_logical[i], true);
            }
            timer.stop();
            validate_sum_degree(sum);
            g_samples.emplace_back("degree_logical_unsorted", num_threads, timer.microseconds());

            g_sum_point_lookups = 0; // reset

            // point lookups, real vertices, sorted
            timer.start();
            sum = 0;
            #pragma omp parallel for num_threads(num_threads) reduction(+:sum) firstprivate(rt, iter_ro)
            for(uint64_t i = 0; i < num_vertices; i++){
                iter_ro.edges(g_vertices_sorted[i], false, [&sum](uint64_t destination, double weight){
                    sum += destination;
                    return false; // stop the iteration
                });
            }
            timer.stop();
            validate_sum_point_lookups(sum);
            g_samples.emplace_back("point_vtx_sorted", num_threads, timer.microseconds());

            // point lookups, real vertices, unsorted
            timer.start();
            sum = 0;
            #pragma omp parallel for num_threads(num_threads) reduction(+:sum) firstprivate(rt, iter_ro)
            for(uint64_t i = 0; i < num_vertices; i++){
                iter_ro.edges(g_vertices_unsorted[i], false, [&sum](uint64_t destination, double weight){
                    sum += destination;
                    return false; // stop the iteration
                });
            }
            timer.stop();
            validate_sum_point_lookups(sum);
            g_samples.emplace_back("point_vtx_unsorted", num_threads, timer.microseconds());

            g_sum_point_lookups = 0; // reset

            // point lookups, logical vertices, sorted
            timer.start();
            sum = 0;
            #pragma omp parallel for num_threads(num_threads) reduction(+:sum) firstprivate(rt, iter_ro)
            for(uint64_t i = 0; i < num_vertices; i++){
                iter_ro.edges(i, true, [&sum](uint64_t destination, double weight){
                    sum += destination;
                    return false; // stop the iteration
                });
            }
            timer.stop();
            validate_sum_point_lookups(sum);
            g_samples.emplace_back("point_logical_sorted", num_threads, timer.microseconds());

            timer.start();
            sum = 0;
            #pragma omp parallel for num_threads(num_threads) reduction(+:sum) firstprivate(rt, iter_ro)
            for(uint64_t i = 0; i < num_vertices; i++){
                iter_ro.edges(g_vertices_logical[i], true, [&sum](uint64_t destination, double weight){
                    sum += destination;
                    return false; // stop the iteration
                });
            }
            timer.stop();
            validate_sum_point_lookups(sum);
            g_samples.emplace_back("point_logical_unsorted", num_threads, timer.microseconds());

            g_sum_scan = 0; // reset

            // scan, real vertices, sorted
            timer.start();
            sum = 0;
            #pragma omp parallel for num_threads(num_threads) reduction(+:sum) firstprivate(rt, iter_ro)
            for(uint64_t i = 0; i < num_vertices; i++){
                iter_ro.edges(g_vertices_sorted[i], false, [&sum](uint64_t destination, double weight){
                    sum += destination;
                    return true; // stop the iteration
                });
            }
            timer.stop();
            validate_sum_scan(sum);
            g_samples.emplace_back("scan_vtx_sorted", num_threads, timer.microseconds());

            // scan, real vertices, unsorted
            timer.start();
            sum = 0;
            #pragma omp parallel for num_threads(num_threads) reduction(+:sum) firstprivate(rt, iter_ro)
            for(uint64_t i = 0; i < num_vertices; i++){
                iter_ro.edges(g_vertices_unsorted[i], false, [&sum](uint64_t destination, double weight){
                    sum += destination;
                    return true; // stop the iteration
                });
            }
            timer.stop();
            validate_sum_scan(sum);
            g_samples.emplace_back("scan_vtx_unsorted", num_threads, timer.microseconds());

            g_sum_scan = 0; // reset

            // scan, logical vertices, sorted
            timer.start();
            sum = 0;
            #pragma omp parallel for num_threads(num_threads) reduction(+:sum) firstprivate(rt, iter_ro)
            for(uint64_t i = 0; i < num_vertices; i++){
                iter_ro.edges(i, true, [&sum](uint64_t destination, double weight){
                    sum += destination;
                    return true; // stop the iteration
                });
            }
            timer.stop();
            validate_sum_scan(sum);
            g_samples.emplace_back("scan_logical_sorted", num_threads, timer.microseconds());


            // scan, logical vertices, unsorted
            timer.start();
            sum = 0;
            #pragma omp parallel for num_threads(num_threads) reduction(+:sum) firstprivate(rt, iter_ro)
            for(uint64_t i = 0; i < num_vertices; i++){
                iter_ro.edges(g_vertices_logical[i], true, [&sum](uint64_t destination, double weight){
                    sum += destination;
                    return true; // stop the iteration
                });
            }
            timer.stop();
            validate_sum_scan(sum);
            g_samples.emplace_back("scan_logical_unsorted", num_threads, timer.microseconds());

        }
    }

    LOG("Experiment completed");
}
#endif

#if defined(HAVE_GRAPHONE)

static void run_graphone(){
    auto* view = create_static_view(get_graphone_graph(), SIMPLE_MASK & PRIVATE_MASK); // global

    const uint64_t num_vertices = g_vertices_sorted.size();
    common::Timer timer;

    for(int r = 0; r < g_num_repetitions; r++){
        LOG("Repetition: " << (r +1) << "/" << g_num_repetitions);
        for(auto num_threads: g_num_threads){
            LOG("    num threads: " << num_threads);

            // vertices, logical identifiers, sorted
            timer.start();
            uint64_t sum = 0;
            #pragma omp parallel for num_threads(num_threads) reduction(+:sum)
            for(uint64_t i = 0; i < num_vertices; i++){
                sum += view->get_degree_out(i);
            }
            timer.stop();
            validate_sum_degree(sum);
            g_samples.emplace_back("degree_logical_sorted", num_threads, timer.microseconds());

            // vertices, logical identifiers, unsorted
            timer.start();
            sum = 0;
            #pragma omp parallel for num_threads(num_threads) reduction(+:sum)
            for(uint64_t i = 0; i < num_vertices; i++){
                sum += view->get_degree_out(g_vertices_logical[i]);
            }
            timer.stop();
            validate_sum_degree(sum);
            g_samples.emplace_back("degree_logical_unsorted", num_threads, timer.microseconds());

            g_sum_point_lookups = 0; // reset

            // point lookups, logical vertices, sorted
            timer.start();
            sum = 0;
            #pragma omp parallel num_threads(num_threads) reduction(+:sum)
            {
                lite_edge_t* neighbours = nullptr;
                uint64_t neighbours_sz = 0;

                #pragma omp for
                for(uint64_t i = 0; i < num_vertices; i++){
                    uint64_t degree = view->get_degree_out(i);
                    if(degree == 0) continue;

                    if(degree > neighbours_sz){
                        neighbours_sz = degree;
                        neighbours = (lite_edge_t*) realloc(neighbours, sizeof(neighbours[0]) * neighbours_sz);
                        if(neighbours == nullptr) throw std::bad_alloc{};
                    }

                    view->get_nebrs_out(i, neighbours);
                    sum += get_sid(neighbours[0]);
                }

                free(neighbours); neighbours = nullptr; neighbours_sz = 0;
            }
            timer.stop();
            validate_sum_point_lookups(sum);
            g_samples.emplace_back("point_logical_sorted", num_threads, timer.microseconds());

            // point lookups, logical, unsorted
            timer.start();
            sum = 0;
            #pragma omp parallel num_threads(num_threads) reduction(+:sum)
            {
                lite_edge_t* neighbours = nullptr;
                uint64_t neighbours_sz = 0;

                #pragma omp for
                for(uint64_t i = 0; i < num_vertices; i++){
                    uint64_t degree = view->get_degree_out(g_vertices_logical[i]);
                    if(degree == 0) continue;

                    if(degree > neighbours_sz){
                        neighbours_sz = degree;
                        neighbours = (lite_edge_t*) realloc(neighbours, sizeof(neighbours[0]) * neighbours_sz);
                        if(neighbours == nullptr) throw std::bad_alloc{};
                    }

                    view->get_nebrs_out(g_vertices_logical[i], neighbours);
                    sum += get_sid(neighbours[0]);
                }

                free(neighbours); neighbours = nullptr; neighbours_sz = 0;
            }
            timer.stop();
            validate_sum_point_lookups(sum);
            g_samples.emplace_back("point_logical_unsorted", num_threads, timer.microseconds());

            // scan, logical vertices, sorted
            timer.start();
            sum = 0;
            #pragma omp parallel num_threads(num_threads) reduction(+:sum)
            {
                lite_edge_t* neighbours = nullptr;
                uint64_t neighbours_sz = 0;

                #pragma omp for
                for(uint64_t i = 0; i < num_vertices; i++){
                    uint64_t degree = view->get_degree_out(i);
                    if(degree == 0) continue;

                    if(degree > neighbours_sz){
                        neighbours_sz = degree;
                        neighbours = (lite_edge_t*) realloc(neighbours, sizeof(neighbours[0]) * neighbours_sz);
                        if(neighbours == nullptr) throw std::bad_alloc{};
                    }

                    view->get_nebrs_out(i, neighbours);

                    for(uint64_t j = 0; j < degree; j++){
                        uint64_t u = get_sid(neighbours[j]);
                        sum += u;
                    }
                }

                free(neighbours); neighbours = nullptr; neighbours_sz = 0;
            }
            timer.stop();
            validate_sum_scan(sum);
            g_samples.emplace_back("scan_logical_sorted", num_threads, timer.microseconds());

            // scan, logical vertices, unsorted
            timer.start();
            sum = 0;
            #pragma omp parallel num_threads(num_threads) reduction(+:sum)
            {
                lite_edge_t* neighbours = nullptr;
                uint64_t neighbours_sz = 0;

                #pragma omp for
                for(uint64_t i = 0; i < num_vertices; i++){
                    uint64_t degree = view->get_degree_out(g_vertices_logical[i]);
                    if(degree == 0) continue;

                    if(degree > neighbours_sz){
                        neighbours_sz = degree;
                        neighbours = (lite_edge_t*) realloc(neighbours, sizeof(neighbours[0]) * neighbours_sz);
                        if(neighbours == nullptr) throw std::bad_alloc{};
                    }

                    view->get_nebrs_out(g_vertices_logical[i], neighbours);

                    for(uint64_t j = 0; j < degree; j++){
                        uint64_t u = get_sid(neighbours[j]);
                        sum += u;
                    }
                }

                free(neighbours); neighbours = nullptr; neighbours_sz = 0;
            }
            timer.stop();
            validate_sum_scan(sum);
            g_samples.emplace_back("scan_logical_unsorted", num_threads, timer.microseconds());
        }
    }

    delete_static_view(view);

    LOG("Experiment completed");
}

#endif

static void validate_sum_degree(uint64_t sum) {
    if(g_sum_degree == 0){
        g_sum_degree = sum;
    } else if (g_sum_degree != sum ){
        cerr << "ERROR: degree mismatch, got: " << sum << ", expected: " << g_sum_degree << "\n";
        throw std::runtime_error("degree mismatch");
    }
}

static void validate_sum_point_lookups(uint64_t sum) {
    if(g_sum_point_lookups == 0){
        g_sum_point_lookups = sum;
    } else if (g_sum_point_lookups != sum ){
        cerr << "ERROR: point lookup sum mismatch, got: " << sum << ", expected: " << g_sum_point_lookups << "\n";
        throw std::runtime_error("point lookup sum mismatch");
    }
}

static void validate_sum_scan(uint64_t sum) {
    if(g_sum_scan == 0){
        g_sum_scan = sum;
    } else if (g_sum_scan != sum ){
        cerr << "ERROR: scan sum mismatch, got: " << sum << ", expected: " << g_sum_scan << "\n";
        throw std::runtime_error("scan mismatch");
    }
}

static void compute_medians(){
    unordered_map</* experiment */ string, unordered_map</* num_threads*/ int, /* completion times */ vector<uint64_t>> > results;

    // collect the results
    for(uint64_t i = 0; i < g_samples.size(); i++){
        auto& sample = g_samples[i];
        results[ sample.m_experiment ][sample.m_num_threads].push_back(sample.m_microsecs);
    }

    // compute the medians
    for(auto& kv: results){
        g_experiments.push_back(kv.first);

        for(auto& kvt: kv.second){
            int num_threads = kvt.first;
            auto& list = kvt.second;

            uint64_t median = 0;
            sort(list.begin(), list.end());
            if(list.size() % 2 == 0){
                median = (list[list.size() /2]  + list[list.size()/2 -1]) /2;
            } else {
                median = list[list.size() /2];
            }


            string key = kv.first;
            key += "@";
            key += num_threads;
            g_medians[key] = median;
        }
    }

    // remove the duplicates from the list of experiments
    sort(g_experiments.begin(), g_experiments.end());
    g_experiments.erase( unique( g_experiments.begin(), g_experiments.end() ), g_experiments.end() );
}

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

    }

    LOG("Loading the graph from " << g_path_graph << " ...");
    auto edges = make_shared<gfe::graph::WeightedEdgeStream> ( g_path_graph );
    edges->permute();

    { // list of vertices
        LOG("Loading the list of vertices...");
        common::Timer timer;
        timer.start();

        // it seems just quicker to read it from the file
        reader::GraphalyticsReader reader ( g_path_graph );
        uint64_t vertex;
        while( reader.read_vertex(vertex) ) { g_vertices_sorted.push_back(vertex); }

        // permute the list
        auto ptr_logical_vertices = make_unique<uint64_t[]>(g_vertices_sorted.size());
        common::permute(ptr_logical_vertices.get(), g_vertices_sorted.size(), /* seed */ 42);
        g_vertices_logical.reserve(g_vertices_sorted.size());
        g_vertices_unsorted.reserve(g_vertices_logical.size());
        for(uint64_t i = 0, end = g_vertices_sorted.size(); i < end; i++){
            g_vertices_logical.emplace_back(ptr_logical_vertices[i]);
            g_vertices_unsorted.emplace_back(g_vertices_sorted[ptr_logical_vertices[i]]);
        }
        sort(g_vertices_sorted.begin(), g_vertices_sorted.end());

        timer.stop();
        LOG("List of vertices loaded in " << timer);
    }


    LOG("Inserting " << edges->num_edges() << " edges into `" << g_library << "' ...");
    gfe::experiment::InsertOnly insert(dynamic_pointer_cast<gfe::library::UpdateInterface>(g_interface), edges, thread::hardware_concurrency(), false);
    insert.execute();
}

static void parse_args(int argc, char* argv[]) {
    static struct option long_options[] = {
        {"graph", required_argument, nullptr, 'G'},
        {"help", no_argument, nullptr, 'h'},
        {"library", required_argument, nullptr, 'l'},
        {"num_threads", required_argument, nullptr, 't'},
        {"repetitions", required_argument, nullptr, 'R'},
        {0, 0, 0, 0} // keep at the end
    };

    int option { 0 };
    int option_index = 0;
    while( (option = getopt_long(argc, argv, "e:G:hl:t:", long_options, &option_index)) != -1 ){
        switch(option){
        case 'e':

            break;
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
            if(library != "teseo" && library != "graphone"){
                cerr << "ERROR: Invalid library: `" << library << "'. Only \"teseo\" and \"graphone\" are supported." << endl;
            }
            g_library = library;
        } break;
        case 'R': {
            g_num_repetitions = stoi(optarg);
            if(g_num_repetitions <= 0){
                cerr << "ERROR: Invalid number of repetitions specified: `" << optarg << "'" << endl;
                exit(EXIT_FAILURE);
            }
        } break;
        case 't': {
            int start = -1;
            int end = 0;
            string str_digit;
            int pos = 0;
            bool stop = false;
            while(!stop){
                if(isdigit(optarg[pos])){
                    str_digit += optarg[pos];
                } else if(optarg[pos] == '-'){
                    if(str_digit.empty() || start > 0){
                        cerr << "ERROR: Invalid format: `" << optarg << "'" << endl;
                        exit(EXIT_FAILURE);
                    }
                    start = stoi(str_digit);
                    end = -1;
                    str_digit.clear();
                } else if(optarg[pos] == ',' || optarg[pos] == '\0'){
                    if(str_digit.empty()){
                        cerr << "ERROR: Invalid format: `" << optarg << "'" << endl;
                        exit(EXIT_FAILURE);
                    }
                    end = stoi(str_digit);

                    if(start < 0){
                        g_num_threads.push_back(end);
                    } else {
                        for(int i = start; i <= end; i++){
                            g_num_threads.push_back(i);
                        }
                    }

                    // next iteration
                    start = -1;
                    end = -1;
                    str_digit.clear();

                    stop = (optarg[pos] == '\0');
                } else if(isspace(optarg[pos])){
                    /* nop, ignore spaces */
                } else {
                    cerr << "ERROR: Invalid format: `" << optarg << "'" << endl;
                    exit(EXIT_FAILURE);
                }

                // next iteration
                pos++;
            }

        } break;
        default:
            assert(0 && "Invalid option");
        }
    }


    if(g_path_graph.empty()){
        cerr << "ERROR: Input graph (-G) not specified\n";
        cerr << string_usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    if(g_num_threads.empty()){
        auto topo = new common::cpu_topology();
        auto t1 = topo->get_threads(/* does not matter */ true, /* smt ? */ false);
        for(uint64_t i = 1; i <= t1.size(); i++){
            g_num_threads.emplace_back(i);
        }
        auto t2 = topo->get_threads(/* does not matter */ true, /* smt ? */ true);
        if(t2.size() > t1.size()){ // it has SMT?
            g_num_threads.emplace_back(t2.size());
        }
        delete topo;
    } else { // remove duplicates
        sort( g_num_threads.begin(), g_num_threads.end() );
        g_num_threads.erase( unique( g_num_threads.begin(), g_num_threads.end() ), g_num_threads.end() );
    }

    if(g_num_repetitions == 0){
        g_num_repetitions = 5;
    }

    if(g_library.empty()){
        g_library = "teseo";
    }
}

static string string_usage(char* program_name) {
    stringstream ss;
    ss << "Usage: " << program_name << " -G <graph> [-t <num_threads>] [-l <library>] [-R <num_repetitions>]\n";
    ss << "Where: \n";
    ss << "  -G <graph> is an .properties file of an undirected graph from the Graphalytics data set\n";
    ss << "  -l <library> is the library to execute. Only \"teseo\" (default) and \"graphone\" are supported\n";
    ss << "  -R <num_repetitions> is the number of repetitions the same micro benchmarks need to be performed\n";
    ss << "  -t <num_threads> follows the page range format, e.g. 1-16,32\n";
    return ss.str();
}

static const int column_sz = 20;

static void print_line(){
    printf("+");
    for(uint64_t i = 0; i < 30; i++){ printf("-"); }
    printf("+");
    for(uint64_t i = 0; i < g_num_threads.size(); i++){
        for(uint64_t i = 0; i < column_sz; i++){ printf("-"); }
        printf("+");
    }
    printf("\n");
}

static void print_header (){
    printf("|%-30s|", " Experiment");
    for(uint64_t i = 0; i < g_num_threads.size(); i++){
        string thread = to_string(g_num_threads[i]);
        thread += " thread(s)";
        printf(" %-*s|", column_sz -1, thread.c_str());
    }
    printf("\n");
}

static int64_t get_median(const std::string& experiment, int num_threads){
    string key = experiment;
    key += "@";
    key += num_threads;

    if(g_medians.count(key) > 0){
        return g_medians[key];
    } else {
        return -1;
    }
}

static string get_scalability(uint64_t t1, uint64_t tn){
    double res = static_cast<double>(t1) / tn;
    char buffer[256];
    sprintf(buffer, "%.2f", res);
    return buffer;
}

static void print_results() {
    print_line();
    print_header();
    print_line();
    for(auto& experiment: g_experiments ){
        printf("| %-29s|", experiment.c_str());

        for(uint64_t i = 0; i < g_num_threads.size(); i++){
            int num_threads = g_num_threads[i];
            string result = to_string(chrono::microseconds{ get_median(experiment, num_threads) } );
            if(i > 0 && g_num_threads[0] == 1){
                result += " (";
                result += get_scalability(get_median(experiment, 1), get_median(experiment, num_threads));
                result += "x)";
            }

            printf("%*s|", column_sz, result.c_str());
        }

        printf("\n");
        print_line();
    }
}

// https://stackoverflow.com/questions/16357999/current-date-and-time-as-string
static string current_date(){
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[256];

    time (&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer,sizeof(buffer),"%d-%m-%Y %H:%M:%S",timeinfo);

    return buffer;
}

static void save_results(const std::string& where) {
    string date = current_date();

    fstream out(where, ios::out);

    out << "{";
    out << "\"version\": 200521, ";
    out << "\"date\": \"" << date << "\", ";
    out << "\"hostname\": \"" << common::hostname() << "\", ";
    out << "\"git_version\": \"" << common::git_last_commit() << "\", ";
    out << "\"samples\": [";
    for(uint64_t i = 0; i < g_samples.size(); i++){
        if(i > 0) out << ", ";
        out << "{";
        out << "\"date\": \"" << date << "\", ";
        out << "\"hostname\": \"" << common::hostname() << "\", ";
        out << "\"library\": \"" << g_library << "\", ";
        out << "\"experiment\": \"" << g_samples[i].m_experiment << "\", ";
        out << "\"num_threads\": \"" << g_samples[i].m_num_threads << "\", ";
        out << "\"microseconds\": \"" << g_samples[i].m_microsecs << "\"";
        out << "}";
    }
    out << "] }";

    out.close();
}


// Time
template <typename D>
static string to_nanoseconds(D duration){
    using namespace std::chrono;
    stringstream result;
    result << (uint64_t) duration_cast<chrono::nanoseconds>(duration).count() << " ns";
    return result.str();
}

template <typename D>
static string to_microseconds(D duration){
    using namespace std::chrono;
    uint64_t time_in_nanosecs = (uint64_t) duration_cast<chrono::nanoseconds>(duration).count();
    uint64_t time_in_microsecs = time_in_nanosecs / 1000;

    stringstream result;
    if(time_in_microsecs >= 3){
        result << time_in_microsecs << " us";
    } else {
        char buffer[128];
        snprintf(buffer, 128, "%.3d", (int) (time_in_nanosecs % 1000));
        result << time_in_microsecs << "." << buffer << " us";
    }

    return result.str();
}

template <typename D>
static string to_milliseconds(D duration){
    using namespace std::chrono;
    uint64_t time_in_microsecs = (uint64_t) duration_cast<chrono::microseconds>(duration).count();
    uint64_t time_in_millisecs = time_in_microsecs / 1000;

    stringstream result;
    if(time_in_microsecs >= 3){
        result << time_in_millisecs << " ms";
    } else {
        char buffer[128];
        snprintf(buffer, 128, "%.3d", (int) (time_in_microsecs % 1000));
        result << time_in_millisecs << "." << buffer << " ms";
    }

    return result.str();
}

template <typename D>
static string to_seconds(D duration){
    using namespace std::chrono;
    uint64_t time_in_millisecs = (uint64_t) duration_cast<chrono::milliseconds>(duration).count();
    uint64_t time_in_seconds = time_in_millisecs / 1000;

    stringstream result;
    char buffer[128];
    snprintf(buffer, 128, "%.3d", (int) (time_in_millisecs % 1000));
    result << time_in_seconds << "." << buffer << " s";

    return result.str();
}

template <typename D>
static string to_minutes(D duration){
    using namespace std::chrono;
    uint64_t seconds = ((uint64_t) duration_cast<chrono::seconds>(duration).count()) % 60ull;
    uint64_t minutes = (uint64_t) duration_cast<chrono::minutes>(duration).count();

    char buffer[128];
    snprintf(buffer, 128, "%" PRIu64 ":%02" PRIu64 " mins", minutes, seconds);
    return string(buffer);
}

template <typename D>
static string to_hours(D duration){
    using namespace std::chrono;
    uint64_t seconds = ((uint64_t) duration_cast<chrono::seconds>(duration).count()) % 60ull;
    uint64_t minutes = (uint64_t) duration_cast<chrono::minutes>(duration).count() % 60ull;
    uint64_t hours = (uint64_t) duration_cast<chrono::hours>(duration).count();

    char buffer[128];
    snprintf(buffer, 128, "%" PRIu64 ":%02" PRIu64 ":%02" PRIu64 " hours", hours, minutes, seconds);
    return string(buffer);
}

template <typename D>
static string to_days(D duration){
    using namespace std::chrono;
    uint64_t seconds = ((uint64_t) duration_cast<chrono::seconds>(duration).count()) % 60ull;
    uint64_t minutes = (uint64_t) duration_cast<chrono::minutes>(duration).count() % 60ull;
    uint64_t hours = (uint64_t) duration_cast<chrono::hours>(duration).count() % 24ull;
    uint64_t days = (uint64_t) duration_cast<chrono::hours>(duration).count() / 24;

    char buffer[128];
    snprintf(buffer, 128, "%" PRIu64 " day(s) and %" PRIu64 ":%02" PRIu64 ":%02" PRIu64 " hours", days, hours, minutes, seconds);
    return string(buffer);
}


template <typename Duration>
static string to_string(Duration d){
    uint64_t time_in_nanosecs = (uint64_t) chrono::duration_cast<chrono::nanoseconds>(d).count();
    if(time_in_nanosecs <= 1000){
        return to_nanoseconds(d);
    } else if(time_in_nanosecs <= (uint64_t) pow(10, 6)){
        return to_microseconds(d);
    } else if(time_in_nanosecs <= (uint64_t) pow(10, 9)) {
        return to_milliseconds(d);
    } else if(time_in_nanosecs <= (uint64_t) pow(10, 9) * 90){ // 90 seconds
        return to_seconds(d);
    } else if(time_in_nanosecs < (uint64_t) pow(10, 9) * 60 * 60){
        return to_minutes(d);
    } else {
        return to_hours(d);
    }
}


