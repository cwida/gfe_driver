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

#include "aging_thread.hpp"

#include <atomic>
#include <random>

#include <gmpxx.h>

#include "common/system.hpp"
#include "experiment/aging.hpp"
#include "graph/edge.hpp"
#include "graph/edge_stream.hpp"
#include "library/interface.hpp"
#include "aging_operation.hpp"
#include "aging_partition.hpp"
#include "async_batch.hpp"
#include "configuration.hpp"

using namespace common;
using namespace std;

/*****************************************************************************
 *                                                                           *
 * Debug                                                                     *
 *                                                                           *
 *****************************************************************************/
#define DEBUG
#define COUT_DEBUG_FORCE(msg) { LOG("[AgingThread::" << __FUNCTION__ << "] [" << concurrency::get_thread_id() << "] " << msg); }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif

/*****************************************************************************
 *                                                                           *
 * AgingThread                                                               *
 *                                                                           *
 *****************************************************************************/

namespace experiment::details {

AgingThread::AgingThread(Aging* instance, const std::vector<AgingPartition>& partitions, int worker_id) : m_instance(instance),
        m_interface(instance->m_interface.get()), m_worker_id(worker_id), m_is_undirected(!m_instance->is_directed()),
        m_partitions(partitions), m_batch(nullptr) {

    if(m_instance->m_batch_size > 0){ // execute updates in batches
        m_batch = new details::AsyncBatch(m_interface, worker_id +1, 4, m_instance->m_batch_size);
    }

    // compute the number of vertices we are responsible to handle
    for(uint64_t part_id = 0, sz = partitions.size(); part_id < sz; part_id++){
        if(m_partitions[part_id].m_start > m_instance->m_max_vertex_id_artificial) break;
        m_num_src_vertices_in_partitions += partitions[part_id].m_length;
    }

    if(is_undirected()){ // compute the max number of edges in each partition
        m_num_edges_in_partition = new mpz_class[m_partitions.size()]; // all mpz instances are set to 0
        mpz_class M { m_instance->m_max_vertex_id_artificial };
        for(uint64_t part_id = 0, sz = m_partitions.size(); part_id < sz; part_id++){
            if(m_partitions[part_id].m_start > m_instance->m_max_vertex_id_artificial) break; // we are not supposed to generate edges for these extra partitions
            mpz_class start { m_partitions[part_id].m_start }; // incl.
            mpz_class end = start + m_partitions[part_id].m_length; // excl.
            mpz_class sum1 = (M - start) * (M - start +1) / 2;
            mpz_class sum2 = (M - end) * (M - end +1) / 2;
            assert(sum1 > sum2);
            m_num_edges_in_partition[part_id] = sum1 - sum2;
        }
    }
}

AgingThread::~AgingThread(){
    delete m_batch; m_batch = nullptr;
    delete[](m_num_edges_in_partition); m_num_edges_in_partition = nullptr;
}

std::future<void> AgingThread::execute(AgingOperation operation){
    auto current_op = m_current_operation;
    if(current_op != AgingOperation::NONE) ERROR("Invalid state: " << (int64_t) operation << ". The worker id " << m_worker_id << " is already busy performing another operation");

    // critical section
    unique_lock<mutex> lock(m_mutex_op);
    m_callback = promise<void>{ };
    auto future = m_callback.get_future();
    m_current_operation = operation;
    lock.unlock();

    m_condvar_op.notify_all();

    return future;
}

void AgingThread::main_thread(){
    while(true){
        unique_lock<mutex> lock(m_mutex_op);
        m_condvar_op.wait(lock, [this](){ return m_current_operation != AgingOperation::NONE; } );

        while(m_current_operation == AgingOperation::NONE) { } // active wait
        if(m_current_operation == AgingOperation::STOP) break; // exit from the loop

        switch(m_current_operation){
        case AgingOperation::START:
            m_interface->on_thread_init(m_worker_id);
            break;
        case AgingOperation::COMPUTE_FINAL_EDGES: {
            m_edges.clear();
            m_final_edges_current_position = 0;
            graph::WeightedEdgeStream* stream = m_instance->m_stream.get();
            for(uint64_t i = 0, end = stream->num_edges(); i < end; i++){
                auto edge = stream->get(i);
                if(is_undirected() && edge.source() > edge.destination()) edge.swap_src_dst(); // ensure src_id < dst_id for undirected graphs
                if(vertex_belongs(edge.source())) m_edges.push_back(edge);
            }
        } break;
        case AgingOperation::EXECUTE_EXPERIMENT: {
            main_experiment();
        } break;
        case AgingOperation::INTERNAL_CLEANUP: {
            m_edges_already_inserted.clear();
        } break;
        case AgingOperation::REMOVE_VERTICES : {
            uint64_t* __restrict vertices = m_instance->m_vertices2remove + m_interval_vertices2remove.m_start;
            for(uint64_t i = 0, sz = m_interval_vertices2remove.m_length; i < sz; i++){
                m_interface->remove_vertex(vertices[i]);
            }
        } break;
        default:
            assert(false && "Invalid operation");
        }

        m_current_operation = AgingOperation::NONE; // move on
        lock.unlock(); // release the lock
        m_callback.set_value();
    }

    assert(m_current_operation == AgingOperation::STOP && "Invalid state, it should still be in the loop");
    m_interface->on_thread_destroy(m_worker_id);
//    m_callback.set_value_at_thread_exit();
    m_callback.set_value();
}

void AgingThread::main_experiment(){
    // constants
    const uint64_t max_number_edges = static_cast<uint64_t>(m_instance->m_ef_edges * m_instance->m_num_edges);
    const int64_t num_total_ops = static_cast<int64_t>(m_instance->m_num_operations_total);
    // heuristics to bump up the probability of inserting a final edge due to multiple threads and deletions
    const double prob_bump = 1.0 * m_instance->m_num_threads;
    const bool report_progress = m_instance->m_report_progress;

    mpz_class max_num_edges = !is_undirected() ? /* ignore */ mpz_class{ 0 } : get_num_edges_in_my_partitions();
    uniform_int_distribution<uint64_t> genrndsrc {0, m_num_src_vertices_in_partitions -1 }; // incl.
    uniform_int_distribution<uint64_t> genrnddst {0, m_instance->m_max_vertex_id_artificial }; // incl
//    uniform_int_distribution<uint64_t> genrndarc {0, max_edge_id > numeric_limits<uint64_t>::max() ?  numeric_limits<uint64_t>::max() : static_cast<uint64_t>(max_edge_id) }; // incl
    gmp_randclass genrndarc(gmp_randinit_default); genrndarc.seed(std::random_device{}());
    uniform_real_distribution<double> genrndweight{0, m_instance->m_max_weight }; // incl.

    int64_t num_ops_done = 0;
    int lastset_coeff = 0;

    if((is_undirected() && max_num_edges > 0) || (!is_undirected() && m_num_src_vertices_in_partitions > 0)){ // edge case
        while( (num_ops_done = m_instance->m_num_operations_performed.fetch_add(m_instance->m_granularity) ) < num_total_ops ){

            // shall we perform a burst of insertions or deletions ?
            if(m_edges2remove.empty() /* There are no edges to remove */ ||
                    m_interface->num_edges() < max_number_edges /* the size of the current graph is no more than (exp_factor)x of the final graph */){

                if(report_progress && static_cast<int>(100.0 * num_ops_done/num_total_ops) > m_instance->m_last_progress_reported){
                    m_instance->m_last_progress_reported = 100.0 * num_ops_done/num_total_ops;
#if defined(DEBUG)
                    LOG("[thread: " << common::concurrency::get_thread_id() << "] "
                            "Progress: " << num_ops_done << "/" << num_total_ops << " (" << 100.0 * num_ops_done/num_total_ops << "%), "
                            "edges final graph: " <<  m_final_edges_current_position << "/" << m_edges.size() << " (" << (100.0 * m_final_edges_current_position/m_edges.size()) << " %)"
                    );
#else // just report the percentages
                    LOG("[thread: " << common::concurrency::get_thread_id() << "] "
                            "Progress: " << 100.0 * num_ops_done/num_total_ops << "%, "
                            "edges final graph: " <<  (100.0 * m_final_edges_current_position/m_edges.size()) << " %"
                    );
#endif
                }

                // insert `m_granularity' edges then
                for(int64_t i = 0, end = m_instance->m_granularity; i < end; i++){
                    double prob_insert_final = prob_bump * static_cast<double>(missing_edges_final()) / (num_total_ops - num_ops_done);
                    if ( m_uniform(m_random) < prob_insert_final){ // insert from the final graph
                        assert(m_final_edges_current_position < m_edges.size());

                        auto edge = m_edges[m_final_edges_current_position];
                        m_final_edges_current_position++;
                        assert((!is_undirected() || edge.source() < edge.destination()) && "Edges in undirected graphs should always be retrieved with the src < dst");

                        // if this edge has been previously inserted remove it
                        auto raw_edge = edge.edge();
                        auto res = m_edges_already_inserted.insert_or_assign(raw_edge, /* final */ true);
                        if(! res.second ){ /* the edge was already present */
                            remove_edge(raw_edge);
                        }

                        insert_edge(edge);

                    } else {
                        // insert a random edge (noise)
                        uint64_t src_id {0}, dst_id {0};
                        double weight = genrndweight(m_random);

                        if(is_undirected()) {
                            mpz_class edge_id = genrndarc.get_z_range(max_num_edges); // generate a random value in [0, max_num_edges)
                            edge_id_2_vertices_id(edge_id, &src_id, &dst_id);
                            assert(src_id < dst_id && "for undirected graphs, any edge should be generated with src < dst");
                        } else {
                            src_id = src_rel2abs(genrndsrc(m_random));
                            dst_id = genrnddst(m_random);
                            if(dst_id == src_id){ // avoid having the same src & dst for an edge
                                dst_id = (dst_id != m_instance->m_max_vertex_id_artificial) ? dst_id +1 : 0;
                            }
                        }
                        assert(src_id <= m_instance->m_max_vertex_id_artificial);
                        assert(dst_id <= m_instance->m_max_vertex_id_artificial);

                        graph::WeightedEdge edge { src_id, dst_id, weight };

                        auto raw_edge = edge.edge();
                        auto res = m_edges_already_inserted.find(raw_edge);
                        bool do_insert = true;
                        bool is_already_registered = false;
                        if(res == m_edges_already_inserted.end()){ // this edge is not already present
                            m_edges_already_inserted[raw_edge] = false;
                        } else if (!res->second){
                            // already present, but it's not final, overwrite its value
                            is_already_registered = true;
                            remove_edge(raw_edge);
                        } else {
                            // already present and final, that is it belongs to the final graph
                            do_insert = false;
                        }

                        if(do_insert){
                            insert_edge(edge);
                            if(!is_already_registered) m_edges2remove.append(raw_edge); /* otherwise already present in the list of edges to remove */
                        }
                    }
                }
            } else {
                // perform a burst of deletions
                for(uint64_t i = 0, end = std::min<uint64_t>(m_instance->m_granularity, m_edges2remove.size()); i < end && !m_edges2remove.empty(); i++){
                    remove_temporary_edge();
                }
            } // end if (burst of insertions or deletions)

            // report how long it took to perform 1x, 2x, ... updates w.r.t. to the size of the final graph
            int aging_coeff = (num_ops_done + m_instance->m_granularity) / m_instance->m_num_edges;
            if(aging_coeff > lastset_coeff){
                if( m_instance->m_last_time_reported.compare_exchange_strong(/* updates lastset_coeff */ lastset_coeff, aging_coeff) ){
                    uint64_t duration = chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - m_instance->m_time_start ).count();
                    m_instance->m_reported_times[aging_coeff -1] = duration;
                }
            }
        } // end while (operation count)

    } // end if (empty list of vertices to cover)

    // insert the missing edges from the final graph
    COUT_DEBUG("Processed edges: " << m_final_edges_current_position << " / " << m_edges.size() << " (" << (100.0 * m_final_edges_current_position/m_edges.size()) << " %)");

    for( ; m_final_edges_current_position < m_edges.size(); m_final_edges_current_position++){
        auto edge = m_edges[m_final_edges_current_position];

        // if this edge has been previously inserted remove it
        auto raw_edge = edge.edge();
        auto res = m_edges_already_inserted.insert_or_assign(raw_edge, /* final */ true);

        if(! res.second ){ /* the edge was already present */
            remove_edge(raw_edge);
        }

        insert_edge(edge);
    }

    // remove all edges that do not belong to the final graph
    while(!m_edges2remove.empty()){
        remove_temporary_edge();
    }

    // done!
    if(m_batch) m_batch->flush(true); // process all pending updates
}

void AgingThread::insert_edge(graph::WeightedEdge edge){
    // be sure that the vertices source & destination are already present
    m_instance->insert_vertex(edge.source());
    m_instance->insert_vertex(edge.destination());

    if(is_undirected() && random01() < 0.5) edge.swap_src_dst(); // noise

    if(m_batch == nullptr){
        // the function returns true if the edge has been inserted. Repeat the loop if it cannot insert the edge as one of
        // the vertices is still being inserted by another thread
        while ( ! m_interface->add_edge(edge) ) { /* nop */ };
    } else {
        m_batch->add_edge(edge);
    }
}

void AgingThread::remove_edge(graph::Edge edge){
    if(is_undirected() && random01() < 0.5) edge.swap_src_dst(); // noise

    if(m_batch == nullptr){
        m_interface->remove_edge(edge);
    } else { // batch updates
        m_batch->remove_edge(edge);
    }
}

double AgingThread::random01() noexcept {
    return m_uniform(m_random);
}

void AgingThread::remove_temporary_edge(){
    assert(!m_edges2remove.empty());

    bool removed = false;
    while(!m_edges2remove.empty() && !removed){
        auto edge = m_edges2remove[0]; m_edges2remove.pop();
        assert(m_edges_already_inserted.find(edge) != m_edges_already_inserted.end() && "This should be in the list of the edges inserted");
        auto it = m_edges_already_inserted.find(edge);
        if(it->second == false /* this is not an edge of the final graph */ ){
            remove_edge(edge);
            m_edges_already_inserted.erase(it);
            removed = true;
        }
    }
}

int64_t AgingThread::missing_edges_final() const {
    assert(m_final_edges_current_position <= m_edges.size());
    return static_cast<int64_t>(m_edges.size() - m_final_edges_current_position);
}

mpz_class AgingThread::get_num_edges_in_my_partitions() const {
    assert(is_undirected() && "Only for undirected graphs");

    mpz_class total;

    for(uint64_t part_id = 0; part_id < m_partitions.size(); part_id++){
        total += m_num_edges_in_partition[part_id];
    }

    return total;
}

void AgingThread::edge_id_2_vertices_id(const mpz_class& edge_id, uint64_t* out_src_id, uint64_t* out_dst_id){
    assert(is_undirected() && "Only for undirected graphs");
    assert(out_src_id != nullptr && out_dst_id != nullptr && "Output argument missing");
    *out_src_id = *out_dst_id = 0; // reset the initial values
    mpz_t tmp_x; mpz_init(tmp_x); // temporary

    bool stop = false;
    mpz_class e = edge_id;
    uint64_t part_id = 0;
    while(!stop && part_id < m_partitions.size()){
        if(e >= m_num_edges_in_partition[part_id]){
            e -= m_num_edges_in_partition[part_id];
            part_id ++;
        } else {
            stop = true;
        }
    }
    assert(stop == true && "The given edge ID is outside the partitions of this worker");

    const mpz_class M = m_instance->m_max_vertex_id_artificial;
    mpz_class v0 = m_partitions[part_id].m_start;
    mpz_class S = (M +1 -v0) * (M -v0) /2;
    // Again we need to solve the inequality S - [(x^2 +x) /2] <= e
    // that is x^2 +x +2(e -s) >= 0 --> [ -1 + sqrt( 1 - 8(e-s) ) ] / 2
//    uint64_t x = ceil( (-1.0 + sqrt(1ll - 8ll*(e-S)) ) / 2.0 );
    mpz_class sqrt_arg = ( S - e ) * 8 +1;
    mpz_sqrt(tmp_x, sqrt_arg.get_mpz_t());
    mpz_sub_ui(tmp_x, tmp_x, 1u);
    mpz_cdiv_q_ui(tmp_x, tmp_x, 2u); // compute the quotient & round up (the `c' in cdiv)
    assert(mpz_class(tmp_x) <= numeric_limits<uint64_t>::max() && "Overflow");

    // convert the result back to uint64_t
    constexpr uint64_t buffer_sz = 256; char buffer[buffer_sz];
    assert(mpz_sizeinbase(tmp_x, 10) <= buffer_sz && "Conversion overflow");
    mpz_get_str(buffer, 10, tmp_x);
    uint64_t x = strtoull(buffer, nullptr, 10);

    // compute the offset from src_id to dst_id
    uint64_t src_id = m_instance->m_max_vertex_id_artificial -x;
    /*bool*/ stop = false;
    do { // we repeat the computation if we obtain a negative offset, it means we went too far when selecting `src_id'
        mpz_class S1 = (M+1 - src_id) * (M - src_id) /2;
        assert(S1 <= S);

        mpz_set(tmp_x, S.get_mpz_t()); // S
        mpz_sub(tmp_x, tmp_x, S1.get_mpz_t()); // S - S1
        mpz_sub(tmp_x, e.get_mpz_t(), tmp_x); // e - [ S - S1 ]

        if(mpz_cmp_si(tmp_x, 0) < 0){ // underflow => negative offset
            src_id --;
        } else {
            stop = true; // done
        }
    } while (!stop);
    assert(mpz_sizeinbase(tmp_x, 10) <= buffer_sz && "Conversion overflow");
    mpz_get_str(buffer, 10, tmp_x);
    uint64_t offset = strtoull(buffer, nullptr, 10);

    uint64_t dst_id = src_id + 1 + offset;

    assert(src_id < dst_id);
    *out_src_id = src_id;
    *out_dst_id = dst_id;

    mpz_clear(tmp_x);
}

uint64_t AgingThread::src_rel2abs(uint64_t relative_vertex_id) const {
    assert(! m_partitions.empty() && "There are no partitions for this worker");
    int64_t count = static_cast<int64_t>(relative_vertex_id);
    int64_t i = 0, sz = m_partitions.size();
    while(i < sz){
        int64_t length = m_partitions[i].m_length; // cast to int64_t
        if(count - length < 0){
            return m_partitions[i].m_start + count;
        } else {
            count -= length;
            i++;
        }
    }

    assert(false && "Invalid vertex ID");
    ERROR("Invalid relative vertex ID: " << relative_vertex_id);
}

bool AgingThread::vertex_belongs(uint64_t vertex_id) const {
    for(auto& p: m_partitions){
        uint64_t start = p.m_start;
        uint64_t end = start + p.m_length;
        if(vertex_id >= start && vertex_id < end){
            return true;
        }
    }

    return false;
}

void AgingThread::set_partition_vertices_to_remove(uint64_t start, uint64_t length){
    COUT_DEBUG("start: " << start << ", length: " << length);
    m_interval_vertices2remove = { start, length };
}

} // namespace


