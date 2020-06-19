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

#include "teseo_lcc.hpp"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <fstream>
#include <future>
#include <iostream>
#include <limits>
#include <mutex>
#include <pthread.h>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include "common/timer.hpp"
#include "utility/timeout_service.hpp"


#include "teseo/context/global_context.hpp"
#include "teseo/memstore/error.hpp"
#include "teseo/memstore/memstore.hpp"
#include "teseo/memstore/scan.hpp"
#include "teseo/transaction/transaction_impl.hpp"
#include "teseo.hpp"

using namespace std;
using namespace teseo;

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG
namespace gfe { extern mutex _log_mutex [[maybe_unused]]; }
#define COUT_CLASS_NAME "?"
#define COUT_DEBUG_FORCE(msg) { std::scoped_lock<std::mutex> lock{::gfe::_log_mutex}; std::cout << "[" << COUT_CLASS_NAME << "::" << __FUNCTION__ << "] " << msg << std::endl; }
#if defined(DEBUG)
    #define COUT_THREAD_SAFE(msg) { std::scoped_lock<std::mutex> lock{::gfe::_log_mutex}; std::cout << msg << std::endl; }
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_THREAD_SAFE(msg)
    #define COUT_DEBUG(msg)
#endif

/*****************************************************************************
 *                                                                           *
 *  Classes                                                                  *
 *                                                                           *
 *****************************************************************************/
namespace gfe::library {

namespace  {

/**
 * Algorithm parameters
 */
static const uint64_t LCC_NUM_WORKERS = thread::hardware_concurrency(); // number of workers / logical threads to use
static constexpr uint64_t LCC_TASK_SIZE = 1ull << 10; // number of vertices processed in each task
//static const uint64_t LCC_NUM_WORKERS = 1; // number of workers / logical threads to use
//static constexpr uint64_t LCC_TASK_SIZE = 2; // number of vertices processed in each task

/**
 * Forward declarations
 */
class LCC_Master;
class LCC_Worker;

class LCC_Master {
    Teseo * const m_teseo; // instance to the database
    Transaction m_transaction; // the transaction providing isolation for the computation
    unordered_map<uint64_t, double>& m_scores; // the final output of the LCC algorithm, a pair <vertex-id, lcc_score>
    unordered_map<uint64_t, atomic<uint64_t>> m_num_triangles; // number of triangles counted so far for the given vertex
    std::atomic<uint64_t> m_next = 0; // counter to select the next task among the workers
    utility::TimeoutService m_timeout; // timer to check whether we are not spending more time than what allocated (1 hour typically)

    // Reserve the space in the hash maps m_score and m_state so that they can be operated concurrently by each thread/worker
    void initialise();

    // Compute the final scores
    void compute_scores();

public:
    // Constructor
    LCC_Master(Teseo* teseo, unordered_map<uint64_t, double>& output, std::chrono::seconds time_budget);

    // Destructor
    ~LCC_Master();

    // Execute the algorithm
    void execute();

    // Select the next window to process, in the form [vertex_start, vertex_end);
    // Return true if a window/task has been fetched, false if there are no more tasks to process
    bool next_task(uint64_t* output_vtx_start /* inclusive */, uint64_t* output_vtx_end /* exclusive */);

    // Retrieve the instance to the database
    Teseo* teseo();

    // Retrieve the number of triangles associated to the given vertex
    atomic<uint64_t>& num_triangles(uint64_t vertex_id);
};

class LCC_Worker {
    LCC_Master* m_master; // handle to the master instance
    Transaction m_transaction;  // the transaction providing isolation for the computation
    Iterator m_iterator; // an iterator attached to the transaction
    thread m_handle; // underlying thread

    // internal state for #process_vertex()
    uint64_t m_n1 =0; // current vertex being processed
    uint64_t m_n2 =0; // current neighbour being visited
    vector<uint64_t> m_neighbours; // neighbours of n1
    uint64_t m_marker =0; // current position in the neighbours vector, to merge shared neighbours
    uint64_t m_num_triangles =0; // current number of triangles found for `n1'

    // Process the given vertex
    void process_vertex(uint64_t vertex_id);

public:
    // Init
    LCC_Worker(LCC_Master* master, const Transaction& transaction);

    // Destructor
    ~LCC_Worker();

    // Main thread
    void execute();

    // Wait for the worker's thread to terminate
    void join();
};

/*****************************************************************************
 *                                                                           *
 *  LCC_Master                                                               *
 *                                                                           *
 *****************************************************************************/
#undef COUT_CLASS_NAME
#define COUT_CLASS_NAME "LCC_Master"

LCC_Master::LCC_Master(Teseo* teseo, unordered_map<uint64_t, double>& output, std::chrono::seconds time_budget)
    : m_teseo(teseo), m_transaction(teseo->start_transaction(/* read only ? */ true)), m_scores(output), m_timeout(time_budget) {
    // initialise();
}

LCC_Master::~LCC_Master() { }

void LCC_Master::initialise(){
    assert(m_scores.size() == 0 && "Already initialised");

    const uint64_t num_vertices = m_transaction.num_vertices();
    m_num_triangles.reserve(num_vertices);

    for(uint64_t i = 0; i < num_vertices; i++){
        uint64_t vertex_id = m_transaction.vertex_id(i);
        m_num_triangles[vertex_id] = 0;
    }
}

void LCC_Master::compute_scores(){
    //common::Timer timer; timer.start();

    m_scores.reserve(m_transaction.num_vertices());
    for(uint64_t i = 0, N = m_transaction.num_vertices(); i < N; i++){
        uint64_t vertex_id = m_transaction.vertex_id(i);
        uint64_t num_triangles = m_num_triangles[vertex_id];
        if(num_triangles == 0){
            m_scores[vertex_id] = 0;
        } else {
            uint64_t degree = m_transaction.degree(i, /* logical ? */ true);
            uint64_t max_num_edges = degree * (degree -1);
            double score = static_cast<double>(num_triangles) / max_num_edges;
            COUT_DEBUG("vertex: " << vertex_id << ", num triangles: " << num_triangles << ", degree: " << degree << ", score: " << score);
            m_scores[vertex_id] = score;
        }
    }

    //timer.stop();
    //COUT_DEBUG_FORCE("compute_scores executed in: " << timer);
}


void LCC_Master::execute() {
    common::Timer timer; timer.start();

    // init the state and the side information for each vertex
    initialise();
    //COUT_DEBUG_FORCE("Initialisation time: " << timer);

    // start the workers
    assert(LCC_NUM_WORKERS >= 1 && "At least one worker should be set");
    vector<LCC_Worker*> workers;
    workers.reserve(LCC_NUM_WORKERS);
    for(uint64_t worker_id = 0; worker_id < LCC_NUM_WORKERS; worker_id++ ){
        workers.push_back(new LCC_Worker(this, m_transaction));
    }

    // wait for the workers to terminate ...
    for(uint64_t worker_id = 0; worker_id < workers.size(); worker_id++ ){
        workers[worker_id]->join();
        delete workers[worker_id]; workers[worker_id] = nullptr;
    }
    if(m_timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }

    compute_scores();
}

bool LCC_Master::next_task(uint64_t* output_vtx_start /* inclusive */, uint64_t* output_vtx_end /* exclusive */) {
    uint64_t logical_start = m_next.fetch_add(LCC_TASK_SIZE); /* return the previous value of m_next */
    uint64_t num_vertices = m_transaction.num_vertices();
    if(logical_start >= num_vertices){
        return false;
    } else {
        uint64_t logical_end = std::min(logical_start + LCC_TASK_SIZE, num_vertices);

        *output_vtx_start = logical_start;
        *output_vtx_end = logical_end;

        return true;
    }
}


Teseo* LCC_Master::teseo() {
    return m_teseo;
}

atomic<uint64_t>& LCC_Master::num_triangles(uint64_t vertex_id) {
    assert(m_num_triangles.count(vertex_id) == 1 && "Vertex not registered");
    return m_num_triangles[vertex_id];
}

/*****************************************************************************
 *                                                                           *
 *  LCC_Worker                                                               *
 *                                                                           *
 *****************************************************************************/

#undef COUT_CLASS_NAME
#define COUT_CLASS_NAME "LCC_Worker"

LCC_Worker::LCC_Worker(LCC_Master* master, const Transaction& transaction) : m_master(master), m_transaction(transaction), m_iterator(m_transaction.iterator()) {
    m_handle = thread { &LCC_Worker::execute, this };
}

LCC_Worker::~LCC_Worker(){ }

void LCC_Worker::execute() {
    COUT_DEBUG("Worker started");

    m_master->teseo()->register_thread();

    uint64_t v_start, v_end;
    while(m_master->next_task(&v_start, &v_end)){
        for(uint64_t v = v_start; v < v_end; v++){
            uint64_t vertex_id = m_transaction.vertex_id(v);
            process_vertex(vertex_id);
        }
    }

    m_iterator.close();
    m_master->teseo()->unregister_thread();

    COUT_DEBUG("Worker terminated");
}

void LCC_Worker::join(){
    m_handle.join();
}

void LCC_Worker::process_vertex(uint64_t n1) {
    COUT_DEBUG("vertex: " << n1);
    m_n1 = n1;
    m_neighbours.clear();
    m_num_triangles = 0; // reset the number of triangles

    m_iterator.edges(n1, false, [this](uint64_t n2, double){
        if(n2 > m_n1) return false; // we're done with n1
        m_n2 = n2;
        m_neighbours.push_back(n2);
        m_marker = 0; // reset the marker

        m_iterator.edges(n2, false, [this](uint64_t n3, double){
            if(n3 > m_n2) return false; // we're done with `n2'
            assert(m_n1 > m_n2 && m_n2 > n3); // we're looking for triangles of the kind c - b - a, with c > b && b > a

            COUT_THREAD_SAFE("  candidate: " << m_n1 << " - " << m_n2 << " - " << n3 << ", marker[" << m_marker << "] = " << m_neighbours[m_marker]);

            if(n3 > m_neighbours[m_marker]){ // merge with m_neighbours
                do {
                    m_marker ++;
                } while(m_marker < m_neighbours.size() && n3 > m_neighbours[m_marker]);
                if(m_marker >= m_neighbours.size()) return false; // there is nothing left to merge
            }

            if(n3 == m_neighbours[m_marker]){ // match !
                COUT_THREAD_SAFE("    match: " << m_n1 << " - " << m_n2 << " - " << n3 << " and " << m_n1 << " - " << n3 << " - " << m_n2);

                m_num_triangles += 2; // we've discovered both n1 - n2 - n3 and n1 - n3 - n2; with n1 > n2 > n3

                // increase the contribution for n2
                m_master->num_triangles(m_n2) += 2;

                // increase the contribution for n3
                m_master->num_triangles(n3) += 2;

                m_marker++;
                if(m_marker >= m_neighbours.size()) return false; // there is nothing left to merge
            }

            return true; // keep scanning

        });

        return true; // keep scanning
    });

    if(m_num_triangles != 0){
        m_master->num_triangles(m_n1) += m_num_triangles;
    }
}

} // anon namespace

/*****************************************************************************
 *                                                                           *
 *  TeseoLCC                                                                 *
 *                                                                           *
 *****************************************************************************/

#undef COUT_CLASS_NAME
#define COUT_CLASS_NAME "TeseoLCC"


TeseoLCC::TeseoLCC(bool is_directed, bool read_only) : TeseoDriver(is_directed, read_only) {

}


void TeseoLCC::lcc(const char* dump2file) {
    unordered_map<uint64_t, double> output;
    Teseo* teseo = reinterpret_cast<Teseo*>(handle_impl());
    teseo->register_thread();
    auto time_budget = m_timeout;

    { // restrict the scope to allow the dtor to clean up
        LCC_Master algorithm ( teseo, output, time_budget );
        algorithm.execute();
    }

    // there is no need to convert from logical IDs to real IDs here, the output already uses real IDs!
    // ...

    // store the results in the given file
    if(dump2file != nullptr){
        COUT_DEBUG("save the results to: " << dump2file);
        fstream handle(dump2file, ios_base::out);
        if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

        for(const auto& keyvaluepair : output){
            handle << keyvaluepair.first << " " << keyvaluepair.second << "\n";
        }

        handle.close();
    }
}


} // namespace

