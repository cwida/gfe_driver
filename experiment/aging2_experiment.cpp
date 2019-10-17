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

#include "aging2_experiment.hpp"

#include <iostream>
#include <mutex>
#include <thread>

#include "configuration.hpp"
#include "common/error.hpp"
#include "details/aging2_master.hpp"
#include "graph/edge_stream.hpp"
#include "library/interface.hpp"

using namespace std;

/*****************************************************************************
 *                                                                           *
 * Debug                                                                     *
 *                                                                           *
 *****************************************************************************/
extern mutex _log_mutex [[maybe_unused]];
//#define DEBUG
#define COUT_DEBUG_FORCE(msg) { scoped_lock<mutex> lock(_log_mutex); cout << "[Aging2Experiment::" << __FUNCTION__ << "] [" << concurrency::get_thread_id() << "] " << msg << endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif

/*****************************************************************************
 *                                                                           *
 * Aging2Experiment                                                          *
 *                                                                           *
 *****************************************************************************/

namespace experiment {

Aging2Experiment::Aging2Experiment(){ }

void Aging2Experiment::set_library(std::shared_ptr<library::UpdateInterface> library) {
    m_library = library;
}

void Aging2Experiment::set_log(const std::string& path){
    m_path_log = path;
}

void Aging2Experiment::set_max_weight(double value){
    if(value <= 0){ INVALID_ARGUMENT("value <= 0: " << value); }
    m_max_weight = value;
}

void Aging2Experiment::set_parallelism_degree(uint64_t num_threads){
    if(num_threads < 1){ INVALID_ARGUMENT("num_threads < 1: " << num_threads); }
    m_num_threads = num_threads;
}

void Aging2Experiment::set_build_frequency(std::chrono::milliseconds millisecs){
    m_build_frequency = millisecs;
}

void Aging2Experiment::set_report_progress(bool value){
    m_report_progress = value;
}

void Aging2Experiment::set_worker_granularity(uint64_t value){
    if(value < 1){ INVALID_ARGUMENT("value < 1: " << value); }
    m_worker_granularity = value;
}

Aging2Result Aging2Experiment::execute(){
    if(m_library.get() == nullptr) ERROR("Library not set. Use #set_library to set it.");
    if(m_path_log.empty()) ERROR("Path to the log file not set. Use #set_log to set it.")

    details::Aging2Master master{*this};
    return master.execute();
}

} // namespace
