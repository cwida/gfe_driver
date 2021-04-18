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

#include "stinger.hpp"

#include <algorithm>
#include <chrono>
#include <cinttypes>
#include <fstream>
#include <limits>
#include <mutex>
#include <queue>
#include <thread>

#include "common/system.hpp"
#include "stinger_core/stinger.h"
extern "C" {
#include "stinger_alg/weakly_connected_components.h"
}
#include "stinger_error.hpp"

using namespace libcuckoo;
using namespace std;

// Macros
#define STINGER reinterpret_cast<struct stinger*>(m_stinger_graph)
#undef CURRENT_ERROR_TYPE
#define CURRENT_ERROR_TYPE ::gfe::library::StingerError

/******************************************************************************
 *                                                                            *
 *  Debug                                                                     *
 *                                                                            *
 *****************************************************************************/
//#define DEBUG
extern mutex _log_mutex [[maybe_unused]];
#define COUT_DEBUG_FORCE(msg) { scoped_lock<mutex> lock(::_log_mutex); std::cout << "[Stinger::" << __FUNCTION__ << "] [" << common::concurrency::get_thread_id() << "] " << msg << std::endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif


/******************************************************************************
 *                                                                            *
 *  Utility functions                                                         *
 *                                                                            *
 ******************************************************************************/
// dump the content to the given file
static void save(vector<pair<uint64_t, int64_t>>& result, const char* dump2file){
    if(dump2file == nullptr) return; // nop
    COUT_DEBUG("save the results to: " << dump2file)

    fstream handle(dump2file, ios_base::out);
    if(!handle.good()) ERROR("Cannot save the result to `" << dump2file << "'");

    for(auto p : result) {
        handle << p.first << " " << p.second << "\n";
    }

    handle.close();
}

/******************************************************************************
 *                                                                            *
 *  Weakly connected components                                               *
 *                                                                            *
 *****************************************************************************/
namespace gfe::library {

void Stinger::wcc(const char* dump2file) {
    // ignore the timeout as we use the impl~ from stinger

    int64_t nv = stinger_max_nv(STINGER);
    COUT_DEBUG("nv: " << nv);
    auto ptr_component_map = make_unique<int64_t[]>(nv); int64_t* component_map = ptr_component_map.get();
    parallel_shiloach_vishkin_components_of_type(STINGER, component_map, /* type, ignore */ 0); // already implemented in Stinger

    // store the final results (if required)
    auto result = to_external_ids(component_map, get_max_num_mappings()); // convert the internal logical IDs into the external IDs
    save(result, dump2file);
}

} // namespace
