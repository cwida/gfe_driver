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
#define DEBUG
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
 *  LCC Implementation                                                       *
 *                                                                           *
 *****************************************************************************/
namespace gfe::library {

namespace  {

/**
 * Forward declarations
 */
class LCC_ControlBlock_LowLevelAPI;
class LCC_ControlBlock_HighLevelAPI;
class LCC_Implementation;
class LCC_VertexState;
class LCC_Watermark;

/**
 * Keep track which is the minimum vertex ID processed by each worker
 */
class LCC_Watermark {
private:
    LCC_Watermark(const LCC_Watermark&) = delete;
    LCC_Watermark& operator=(const LCC_Watermark&) = delete;

    const uint64_t m_num_workers;
    uint64_t* m_watermarks;

public:
    LCC_Watermark(uint64_t num_workers) : m_num_workers(num_workers){
        m_watermarks = new uint64_t[m_num_workers]();
    }

    ~LCC_Watermark(){
        delete[] m_watermarks;
    }

    void set_watermark(uint64_t worker_id, uint64_t min_vertex_id){
        assert(worker_id < m_num_workers && "Overflow");
        m_watermarks[worker_id] = min_vertex_id;
    }

    uint64_t get_watermark() const {
        uint64_t watermark = numeric_limits<uint64_t>::max();
        for(uint64_t i = 0; i < m_num_workers; i++){
            uint64_t value = m_watermarks[i];
            watermark = min(value, watermark);
        }

        return watermark;
    }
};


/**
 * The number of triangles discovered by visiting a previous vertex (vertex_id)
 */
struct PartialCount {
    uint64_t m_vertex_id; // the source of the triangles
    uint64_t m_contribution; // how many triangles have been discovered starting from source

    PartialCount(uint64_t vertex_id, uint64_t contribution = 0) : m_vertex_id(vertex_id), m_contribution(contribution) { }
};

// Print to the stream the content of the PartialCount instance, for debugging purposes
[[maybe_unused]] static
ostream& operator<<(ostream& out, const PartialCount& pc){
    out << "{vertex_id: " << pc.m_vertex_id << ", contribution: " << pc.m_contribution << "}";
    return out;
}

/**
 * Side information associated to a each vertex to aid the computation
 */
class LCC_VertexState {
    pthread_spinlock_t m_latch; // provide thread safety
    bool m_closed; // check whether other workers are allowed to push new vertices into this vertex
    vector<PartialCount> m_counts; // triangles already computed for this vertex

public:
    // Init
    LCC_VertexState(){
        pthread_spin_init(&m_latch, /* share among multiple processes ? */ false);
        m_closed = false;
    }

    // Destructor
    ~LCC_VertexState(){
        pthread_spin_destroy(&m_latch);
    }

    void insert(PartialCount count){
        // test, lock and test again
        if(!m_closed){
            pthread_spin_lock(&m_latch);
            if(!m_closed){
                m_counts.push_back(count);
            }

            pthread_spin_unlock(&m_latch);
        }
    }

    vector<PartialCount> done(uint64_t watermark){
        vector<PartialCount> result;
        pthread_spin_lock(&m_latch);
        m_closed = true;
        pthread_spin_unlock(&m_latch);

        for(uint64_t i = 0; i < m_counts.size(); i ++){
            if(m_counts[i].m_vertex_id <= watermark){
                result.push_back(m_counts[i]);
            }
        }

        sort(result.begin(), result.end(), [](const PartialCount& p1, const PartialCount& p2){
            return p1.m_vertex_id < p2.m_vertex_id;
        });

        // release the associated memory
        //m_counts.clear();
        vector<PartialCount>().swap(m_counts);

        return result;
    }
};


#undef COUT_CLASS_NAME
#define COUT_CLASS_NAME "LCC_Implementation"

/**
 * LCC kernel implementation
 */
class LCC_Implementation {
    friend class LCC_ControlBlock_HighLevelAPI; // user API
    friend class LCC_ControlBlock_LowLevelAPI; // internal API

private:
    const bool m_is_low_level_api; // whether to use the low level or the user API to implement the kernel
    const uint64_t m_num_workers; // total number of threads used for the computation
    Teseo * const m_teseo; // instance to the database
    Transaction m_transaction; // the transaction providing isolation for the computation
    const uint64_t m_num_vertices; // total number of vertices in the graph
    unordered_map<uint64_t, double>& m_scores; // the final output of the LCC algorithm, a pair <vertex-id, lcc_score>
    unordered_map<uint64_t, LCC_VertexState> m_state; // side information mantained for each vertex, during the computation
    LCC_Watermark m_watermarks; // keep track of the progress of each worker, the recycle previous scores
    constexpr static uint64_t m_step_size = 1024; // the granularity of each task performed by the workers
    std::atomic<uint64_t> m_next = 0; // counter to select the next task among the workers
    utility::TimeoutService m_timeout; // timer to check whether we are not spending more time than what allocated (1 hour typically)

    // Ensure that the current thread is registered in the database
    static Teseo* register_thread(Teseo* teseo){
        teseo->register_thread();
        return teseo;
    }

    // Select the next window to process, in the form [vertex_start, vertex_end);
    // Return true if a window/task has been fetched, false if there are no more tasks to process
    bool next_task(uint64_t* output_vtx_start /* inclusive */, uint64_t* output_vtx_end /* exclusive */, bool logical_vertices){
        assert(output_vtx_start != nullptr && output_vtx_end != nullptr && "Null argument");
        uint64_t logical_start = m_next.fetch_add(m_step_size); /* return the previous value of m_next */
        if(logical_start >= m_num_vertices) {
            return false;
        } else {

            uint64_t logical_end = logical_start + m_step_size;


            if(logical_vertices){
                *output_vtx_start = logical_start;
                *output_vtx_end = std::min(m_num_vertices, logical_end);
            } else {

                *output_vtx_start = m_transaction.vertex_id(logical_start);

                if(logical_end < m_num_vertices){
                    *output_vtx_end = m_transaction.vertex_id(logical_end);
                } else {
                    *output_vtx_end = numeric_limits<uint64_t>::max();
                }

            }

            return true;
        }
    }

    // Process the task for the vertices in [window_start, window_end)
    // Need to be defined afterwards due to the dependency to LCC_ControlBlock_LowLevelAPI
    void process_window_lla(int worker_id, uint64_t window_start, uint64_t window_end); // real vertices
    void process_window_hla(int worker_id, uint64_t window_start, uint64_t window_end); // logical interval

    // Main loop for each worker. Fetch & process one task at the time
    void thread_main(int worker_id){
        // no need for thread pinning

        m_teseo->register_thread();

        uint64_t window_start = 0, window_end = 0;
        if(m_is_low_level_api){
            while(next_task(&window_start, &window_end, /* logical ? */ false) && /* there is still time to process */ !m_timeout.is_timeout()){
                process_window_lla(worker_id, window_start, window_end);
            }
        } else { // user API
            while(next_task(&window_start, &window_end, /* logical ? */ true) && /* there is still time to process */ !m_timeout.is_timeout()){
                process_window_hla(worker_id, window_start, window_end);
            }
        }

        m_teseo->unregister_thread();
    }

    // Reserve the space in the hash maps m_score and m_state so that they can be operated concurrently by each thread/worker
    void initialise(){
        assert(m_scores.size() == 0 && "Already initialised");

        const uint64_t num_vertices = m_transaction.num_vertices();
        m_scores.reserve(num_vertices);
        m_state.reserve(num_vertices);

        for(uint64_t i = 0; i < num_vertices; i++){
            uint64_t vertex_id = m_transaction.vertex_id(i);
            m_scores[vertex_id] = 0.0;
            m_state[vertex_id];
        }
    }

public:
    // Constructor, init the private fields
    LCC_Implementation(bool is_low_level_api, Teseo* teseo, unordered_map<uint64_t, double>& output, std::chrono::seconds time_budget) :
        m_is_low_level_api(is_low_level_api),
        m_num_workers(thread::hardware_concurrency()),
        m_teseo(teseo),
        m_transaction(register_thread(teseo)->start_transaction(/* read only ?*/ true)),
        m_num_vertices(m_transaction.num_vertices()),
        m_scores(output),
        m_watermarks(m_num_workers),
        m_timeout(time_budget) {
    }

    // Destructor
    ~LCC_Implementation(){
        m_transaction.commit(); // rollback or commit, it doesn't really matter for read only transactions
        m_teseo->unregister_thread();
    }

    // Execute the LCC algorithm. Store its output into the hash table _output_ (m_score) provided in the instance
    // of this class
    void execute(){
        common::Timer timer; timer.start();

        // init the state and the side information for each vertex
        initialise();

        // fire the computation
        vector<thread> workers;
        for(int worker_id = 0, N = m_num_workers; worker_id < N; worker_id++){
            workers.emplace_back(&LCC_Implementation::thread_main, this, worker_id);
        }

        // and wait for its completion
        for(auto& t: workers) t.join();

        if(m_timeout.is_timeout()){ RAISE_EXCEPTION(TimeoutError, "Timeout occurred after " << timer); }
    }

};

#undef COUT_CLASS_NAME
#define COUT_CLASS_NAME "LCC_ControlBlock_LowLevelAPI"

/**
 * Control block, process one vertex or one at edge at the time, as their retrieved from the storage
 */
class LCC_ControlBlock_LowLevelAPI {
    LCC_ControlBlock_LowLevelAPI(const LCC_ControlBlock_LowLevelAPI&) = delete;
    LCC_ControlBlock_LowLevelAPI& operator=(const LCC_ControlBlock_LowLevelAPI&) = delete;

    const int m_worker_id; // the ID of the current worker, associated to the underlying thread
    const uint64_t m_last_vertex_id; // the last vertex of the range (exclusive) we need to process in this task
    LCC_Implementation* m_instance; // owner instance of this control block
    vector<PartialCount> m_already_processed; // list of neighbours already processed (vertex id <= watermark)
    vector<PartialCount> m_neighbours; // list of neighbours of the current source (vertex_id > watermark)
    uint64_t m_watermark; // current threshold for the vertices already processed by this or other workers
    uint64_t m_n1; // current source vertex being processed, external user ID
    uint64_t m_num_triangles; // number of triangles counted for n1
    uint64_t m_degree; // number of neighbours/edges attached to n1, counted so far
    Iterator m_iterator; // an iterator attached to the transaction

public:
    // Init
    LCC_ControlBlock_LowLevelAPI(LCC_Implementation* lcc_impl, uint64_t last_vertex_id, int worker_id) :
        m_worker_id(worker_id), m_last_vertex_id(last_vertex_id),
        m_instance(lcc_impl),
        m_watermark(lcc_impl->m_watermarks.get_watermark()),
        m_n1(numeric_limits<uint64_t>::max()), m_num_triangles (0), m_degree(0),
        m_iterator(lcc_impl->m_transaction.iterator()) { }

    // Destructor
    ~LCC_ControlBlock_LowLevelAPI(){
        flush();
    }

    /**
     * Compute the LCC score for the current vertex (m_n1) and propagate the number of computed triangles
     * to the other vertices.
     */
    void flush(){
        if(m_n1 == numeric_limits<uint64_t>::max()) return; // nothing to save

        // propagate the computed triangles to the other vertices
        uint64_t i = 0, sz = m_neighbours.size();
        while(i < sz && m_neighbours[i].m_vertex_id < m_n1) i++; // we missed the watermark
        while(i < sz){
            auto vertex_id = m_neighbours[i].m_vertex_id;
            auto contribution = m_neighbours[i].m_contribution;
            if(contribution > 0){
                assert(m_instance->m_state.count(vertex_id) == 1 && "Entry not initialised");
                m_instance->m_state[vertex_id].insert(PartialCount{ m_n1, contribution });
            }
            i++; // next iteration
        }

        // compute the score for the vertex `m_n1'
        if(m_degree >= 2){ // for degree < 2, the score is 0 by definition
            uint64_t max_num_edges = m_degree * (m_degree -1);
            double score = static_cast<double>(m_num_triangles) / max_num_edges;
            COUT_DEBUG("vertex: " << m_n1 << ", num triangles: " << m_num_triangles << ", degree: " << m_degree << ", score: " << score);
            m_instance->m_scores[m_n1] = score;
        }

        // update the watermark for the current worker
        m_instance->m_watermarks.set_watermark(m_worker_id, m_n1);
    }

    /**
     * Check whether we have exhausted the current window ?
     */
    bool done() {
        if(m_n1 >= m_last_vertex_id){
            m_n1 = numeric_limits<uint64_t>::max();
            return true;
        } else {
            return false;
        }
    }

    /**
     * Process the given vertex, as retrieved from the storage
     */
    bool process_vertex(uint64_t vertex_id){
        COUT_THREAD_SAFE("-- vertex: " << vertex_id);

        flush(); // save the score for the previous vertex
        m_n1 = vertex_id; // update the current vertex being processed
        if(done()) return false; // are we done with the current window, for this worker ?

        // clean the state to process m_n1
        m_neighbours.clear();
        assert(m_instance->m_state.count(m_n1) == 1 && "Entry not initialised");
        m_already_processed = m_instance->m_state[m_n1].done(m_watermark);
        m_degree = 0;
        m_num_triangles = 0;
        for(uint64_t i = 0, sz = m_already_processed.size(); i < sz; i++){
            m_num_triangles += m_already_processed[i].m_contribution;
        }

#if defined(DEBUG)
        {
            std::scoped_lock<std::mutex> lock{::gfe::_log_mutex};
            cout << "  watermark: " << m_watermark << "\n";
            cout << "  already processed: [";
            for(uint64_t i = 0; i < m_already_processed.size(); i++){
                if(i > 0) cout << ", ";
                cout << m_already_processed[i];
            }
            cout << "]\n";
            cout << "  num triangles (sum of already processed): " << m_num_triangles << "\n";
        }
#endif

        return true; // keep going
    }

    /**
     * Process the given edge, that is m_n1 -> n2
     */
    bool process_edge(uint64_t n2){
        // we're trying to build the triangles `n1 - n2 - ?'
        COUT_THREAD_SAFE("    [" << m_degree << "] edge: " << m_n1 << " -> " << n2 << ", watermark: " << m_watermark);
        m_degree++;

        if(n2 <= m_watermark) return true; // we have already processed this vertex
        m_neighbours.emplace_back(n2, /* num of triangles via n2 */ 0); // keep track of our neighbours
        uint64_t marker_wtm = 0; // cursor for the array m_already_processed
        uint64_t marker_nb = 0; // cursor for the array m_neighbours

#if defined(DEBUG)
        {
            std::scoped_lock<std::mutex> lock{::gfe::_log_mutex};
            cout << "        neighbours [";
            for(uint64_t i = 0; i < m_neighbours.size(); i++){
                if(i > 0) cout << ", ";
                cout << i << ": " << m_neighbours[i];
            }
            cout << "]\n";
        }
#endif
        // iterate over the neighbours of n2
        m_iterator.edges(n2, false, [this, n2, &marker_wtm, &marker_nb](uint64_t n3, double){
            if (n3 <= m_watermark) { // merge with the array m_already_processed
                COUT_THREAD_SAFE("        candidate (wtm): " << m_n1 << " - " << n2 << " - " << n3);

                while(marker_wtm < m_already_processed.size() && n3 > m_already_processed[marker_wtm].m_vertex_id){ marker_wtm++; }
                // m_already_processed contains all neighbours of n1 less than `watermark'
                if(marker_wtm < m_already_processed.size() && n3 == m_already_processed[marker_wtm].m_vertex_id){
                    COUT_THREAD_SAFE("          match (wtm): " << marker_wtm);
                    m_num_triangles ++;

                    // increase the contribution for n2
                    assert(m_neighbours.back().m_vertex_id == n2);
                    m_neighbours.back().m_contribution ++;

                    marker_wtm++; // next iteration
                }

                return true; // keep scanning
            } else if(n3 >= n2){
                assert(n3 >= m_watermark); // as n2 >= m_watermark
                return false; // stop the iterator
            } else { // n3 < n2
                COUT_THREAD_SAFE("        candidate (nb): " << m_n1 << " - " << n2 << " - " << n3 << ", marker_nb: " << marker_nb << " = " << m_neighbours[marker_nb]);

                if(n3 > m_neighbours[marker_nb].m_vertex_id){ // merge with m_neighbours
                    do {
                        marker_nb ++;
                    } while(marker_nb < m_neighbours.size() && n3 > m_neighbours[marker_nb].m_vertex_id);
                    if(marker_nb >= m_neighbours.size()) return false; // there is nothing left to merge
                }
                if(n3 == m_neighbours[marker_nb].m_vertex_id){
                    COUT_THREAD_SAFE("          match (nb): " << m_n1 << " - " << n2 << " - " << n3 << " and " << m_n1 << " - " << n3 << " - " << n2);

                    m_num_triangles += 2; // we've discovered both n1 - n2 - n3 and n1 - n3 - n2; n2 > n3

                    // increase the contribution for n2
                    assert(m_neighbours.back().m_vertex_id == n2);
                    m_neighbours.back().m_contribution ++;

                    // increase the contribution for n3
                    m_neighbours[marker_nb].m_contribution ++;

                    marker_nb++;
                    if(marker_nb >= m_neighbours.size()) return false; // there is nothing left to merge
                }

                return true; // keep scanning
            }
        }); // n2.iterator

        return true; // keep going
    }

    /**
     * Process the given vertex/edge. Source and destination are external user IDs
     */
    bool process(bool is_vertex, uint64_t source, uint64_t destination){
        if(is_vertex){
            return process_vertex(source);
        } else { // this is an edge
            assert(source == m_n1 && "We haven`t seen the source vertex");
            return process_edge(destination);
        }
    }

    /**
     * Callback for the scan
     */
    bool operator()(uint64_t source_internal_id, uint64_t destination_internal_id, double){
        return process(/* is vertex ? */ destination_internal_id == 0, /* I2E */ source_internal_id -1, /* I2E */ destination_internal_id -1);
    }
};

#undef COUT_CLASS_NAME
#define COUT_CLASS_NAME "LCC_Implementation"

void LCC_Implementation::process_window_lla(int worker_id, uint64_t window_start, uint64_t window_end){
    assert(m_is_low_level_api == true && "To be used for the low level API");

    // window_start and window_end are already real user IDs
    COUT_DEBUG("worker_id: " << worker_id << ", window: [" << window_start << ", " << window_end << ")");

    auto transaction = m_transaction;

    auto tx_impl = reinterpret_cast<transaction::TransactionImpl*>(transaction.handle_impl());
    auto memstore = teseo::context::global_context()->memstore();
    LCC_ControlBlock_LowLevelAPI cb(this, window_end, worker_id);

    memstore->scan(tx_impl, /* E2I */ window_start +1, 0, cb);
}

#undef COUT_CLASS_NAME
#define COUT_CLASS_NAME "LCC_ControlBlock_HighLevelAPI"

class LCC_ControlBlock_HighLevelAPI {
    LCC_ControlBlock_HighLevelAPI(const LCC_ControlBlock_HighLevelAPI&) = delete;
    LCC_ControlBlock_HighLevelAPI& operator=(const LCC_ControlBlock_HighLevelAPI&) = delete;

    const int m_worker_id; // the ID of the current worker, associated to the underlying thread
    LCC_Implementation* const m_instance; // owner of this control block
    Transaction m_transaction; // the transaction providing isolation to the computation
    Iterator m_iterator; // an iterator attached to the transaction
    uint64_t m_n1 = 0; // current vertex being examined
    uint64_t m_watermark = 0; // current watermark,
    uint64_t m_num_triangles = 0; // number of triangles found for the current vertex
    vector<PartialCount> m_already_processed; // list of neighbours already processed (vertex id <= watermark)
    vector<PartialCount> m_neighbours; // list of neighbours of the current source (vertex_id > watermark)


    void process_vertex(uint64_t n1){ // real vertex ID
        COUT_DEBUG("vertex: " << n1);

        const uint64_t degree = m_transaction.degree(n1, /* logical ? */ false);

        if(degree >= 2) {
            m_watermark = m_instance->m_watermarks.get_watermark();
            m_n1 = n1;
            m_already_processed =  m_instance->m_state[n1].done(m_watermark);
            m_neighbours.clear();
            m_num_triangles = 0;
            for(uint64_t i = 0, sz = m_already_processed.size(); i < sz; i++){
                m_num_triangles += m_already_processed[i].m_contribution;
            }

            m_iterator.edges(n1, false, [this](uint64_t n2, double){
                if(n2 <= m_watermark) return true; // we have already processed this vertex
                m_neighbours.emplace_back(n2, /* num of triangles via n2 */ 0); // keep track of our neighbours
                uint64_t marker_wtm = 0; // cursor for the array m_already_processed
                uint64_t marker_nb = 0; // cursor for the array m_neighbours

                // we are looking for triangles in the route m_n1 - n2 - ?

                // iterate over the neighbours of n2
                m_iterator.edges(n2, false, [this, n2, &marker_wtm, &marker_nb](uint64_t n3, double){
                    if (n3 <= m_watermark) { // merge with the array m_already_processed
                        COUT_THREAD_SAFE("        candidate (wtm): " << m_n1 << " - " << n2 << " - " << n3);

                        while(marker_wtm < m_already_processed.size() && n3 > m_already_processed[marker_wtm].m_vertex_id){ marker_wtm++; }
                        // m_already_processed contains all neighbours of n1 less than `watermark'
                        if(marker_wtm < m_already_processed.size() && n3 == m_already_processed[marker_wtm].m_vertex_id){
                            COUT_THREAD_SAFE("          match (wtm): " << marker_wtm);
                            m_num_triangles ++;

                            // increase the contribution for n2
                            assert(m_neighbours.back().m_vertex_id == n2);
                            m_neighbours.back().m_contribution ++;

                            marker_wtm++; // next iteration
                        }

                        return true; // keep scanning
                    } else if(n3 >= n2){
                        assert(n3 >= m_watermark); // as n2 >= m_watermark
                        return false; // stop the iterator
                    } else { // n3 < n2
                        COUT_THREAD_SAFE("        candidate (nb): " << m_n1 << " - " << n2 << " - " << n3 << ", marker_nb: " << marker_nb << " = " << m_neighbours[marker_nb]);

                        if(n3 > m_neighbours[marker_nb].m_vertex_id){ // merge with m_neighbours
                            do {
                                marker_nb ++;
                            } while(marker_nb < m_neighbours.size() && n3 > m_neighbours[marker_nb].m_vertex_id);
                            if(marker_nb >= m_neighbours.size()) return false; // there is nothing left to merge
                        }
                        if(n3 == m_neighbours[marker_nb].m_vertex_id){
                            COUT_THREAD_SAFE("          match (nb): " << m_n1 << " - " << n2 << " - " << n3 << " and " << m_n1 << " - " << n3 << " - " << n2);

                            m_num_triangles += 2; // we've discovered both n1 - n2 - n3 and n1 - n3 - n2; n2 > n3

                            // increase the contribution for n2
                            assert(m_neighbours.back().m_vertex_id == n2);
                            m_neighbours.back().m_contribution ++;

                            // increase the contribution for n3
                            m_neighbours[marker_nb].m_contribution ++;

                            marker_nb++;
                            if(marker_nb >= m_neighbours.size()) return false; // there is nothing left to merge
                        }

                        return true; // keep scanning
                    }
                }); // n2.iterator


                return true;
            });

            // propagate the computed triangles to the other vertices
            uint64_t i = 0, sz = m_neighbours.size();
            while(i < sz && m_neighbours[i].m_vertex_id < n1) i++; // we missed the watermark
            while(i < sz){
                auto vertex_id = m_neighbours[i].m_vertex_id;
                auto contribution = m_neighbours[i].m_contribution;
                if(contribution > 0){
                    assert(m_instance->m_state.count(vertex_id) == 1 && "Entry not initialised");
                    m_instance->m_state[vertex_id].insert(PartialCount{ n1, contribution });
                }
                i++; // next iteration
            }

            // compute the score for the vertex `n1'
            uint64_t max_num_edges = degree * (degree -1);
            double score = static_cast<double>(m_num_triangles) / max_num_edges;
            COUT_DEBUG("vertex: " << n1 << ", num triangles: " << m_num_triangles << ", degree: " << degree << ", score: " << score);
            m_instance->m_scores[n1] = score;
        }

        // update the watermark for the current worker
        m_instance->m_watermarks.set_watermark(m_worker_id, n1);
    }

public:
    LCC_ControlBlock_HighLevelAPI(LCC_Implementation* instance, int worker_id, Transaction txn) :
        m_worker_id(worker_id), m_instance(instance), m_transaction(txn), m_iterator(m_transaction.iterator()){

    }

    void execute(uint64_t window_start /* incl */, uint64_t window_end /* excl */){
        for(uint64_t i = window_start; i < window_end; i++){
            uint64_t vertex_id = m_transaction.vertex_id(i);
            process_vertex(vertex_id);
        }
    }
};

#undef COUT_CLASS_NAME
#define COUT_CLASS_NAME "LCC_Implementation"

void LCC_Implementation::process_window_hla(int worker_id, uint64_t window_start, uint64_t window_end){
    assert(m_is_low_level_api == false && "To be used with the user API");

    // window_start and window_end are logical IDs
    COUT_DEBUG("worker_id: " << worker_id << ", window: [" << window_start << ", " << window_end << ")");

    LCC_ControlBlock_HighLevelAPI cb(this, worker_id, m_transaction);
    cb.execute(window_start, window_end);
}


} // anon namespace

/*****************************************************************************
 *                                                                           *
 *  TeseoLCC                                                                 *
 *                                                                           *
 *****************************************************************************/

#undef COUT_CLASS_NAME
#define COUT_CLASS_NAME "TeseoLCC"


TeseoLCC::TeseoLCC(bool is_directed, bool low_level_api) : TeseoDriver(is_directed, /* read only ? */ true), m_low_level_api(low_level_api) {

}


void TeseoLCC::lcc(const char* dump2file) {
    unordered_map<uint64_t, double> output;
    Teseo* teseo = reinterpret_cast<Teseo*>(handle_impl());
    auto time_budget = m_timeout;

    { // restrict the scope to allow the dtor to clean up
        LCC_Implementation implementation(m_low_level_api, teseo, output, time_budget);
        implementation.execute();
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

