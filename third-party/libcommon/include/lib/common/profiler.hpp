/**
 * This file is part of libcommon.
 *
 * libcommon is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, orF
 * (at your option) any later version.
 *
 * libcommon is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with libcommon.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef COMMON_PROFILER_HPP
#define COMMON_PROFILER_HPP

#include <cinttypes>
#include <cstddef>
#include <ostream>

#include "error.hpp"
#include "details/profiler.hpp"
#include "database.hpp"

/*****************************************************************************
 *                                                                           *
 *   Exceptions                                                              *
 *                                                                           *
 *****************************************************************************/
namespace common {
DEFINE_EXCEPTION(ProfilerError);
}

/*****************************************************************************
 *                                                                           *
 *   Cache faults                                                            *
 *                                                                           *
 *****************************************************************************/
namespace common {
/**
 * Data recorded by the CachesProfiler
 */
struct CachesSnapshot {
    uint64_t m_cache_l1_misses = 0; // number of misses in the L1
    uint64_t m_cache_llc_misses = 0; // number of misses in the LLC (=L3 assumed)
    uint64_t m_cache_tlb_misses = 0; // number of misses in the TLB (I think from LLC, assumes page-walk)

    void operator+=(CachesSnapshot snapshot);
    Database::BaseRecord data_record() const;
};

/**
 * Record the amount of L1, LLC and TLB faults hit during the execution.
 * Usage:
 *      CachesProfiler profiler;
 *      profiler.start();
 *      ... computation ...
 *      profiler.stop();
 *      results = profiler.snapshot();
 */
class CachesProfiler : public details::GenericProfiler {
    CachesSnapshot m_current_snapshot;

public:
    /**
     * Initialise the profiler
     */
    CachesProfiler();

    /**
     * Destructor
     */
    ~CachesProfiler();


    /**
     * Start recording
     */
    void start();

    /**
     * Retrieve the data associated to this profiler
     */
    CachesSnapshot snapshot();

    /**
     * Stop the recording
     */
    CachesSnapshot stop();

    /**
     * Retrieve a data record ready to be stored in the database
     */
    Database::BaseRecord data_record();
};

std::ostream& operator<<(std::ostream& out, const CachesSnapshot& snapshot);
CachesSnapshot operator+(const CachesSnapshot& s1, const CachesSnapshot& s2);

/*****************************************************************************
 *                                                                           *
 *   Branch mispredictions                                                   *
 *                                                                           *
 *****************************************************************************/

/**
 * Data recorded by the BranchMispredictionsProfiler
 */
struct BranchMispredictionsSnapshot{
    uint64_t m_conditional_branches = 0; // total number of conditional branch instructions
    uint64_t m_branch_mispredictions = 0; // total number of branch mispredictions
    uint64_t m_cache_l1_misses = 0; // number of cache misses in the L1
    uint64_t m_cache_llc_misses = 0; // number of cache misses in the LLC (=L3 assumed)

    void operator+=(BranchMispredictionsSnapshot snapshot);
    Database::BaseRecord data_record() const;
};

/**
 * Record the amount of branch mispredictions AND L1 and LLC faults during the executions
 */
struct BranchMispredictionsProfiler : public details::GenericProfiler {
    BranchMispredictionsSnapshot m_current_snapshot;

public:
    /**
 * Initialise the profiler
 */
    BranchMispredictionsProfiler();

    /**
     * Destructor
     */
    ~BranchMispredictionsProfiler();


    /**
     * Start recording
     */
    void start();

    /**
     * Retrieve the data associated to this profiler
     */
    BranchMispredictionsSnapshot snapshot();

    /**
     * Stop the recording
     */
    BranchMispredictionsSnapshot stop();

    /**
     * Retrieve a data record ready to be stored in the database
     */
    Database::BaseRecord data_record();
};

std::ostream& operator<<(std::ostream& out, const BranchMispredictionsSnapshot& snapshot);


/*****************************************************************************
 *                                                                           *
 *   SoftwareEventsProfiler                                                  *
 *                                                                           *
 *****************************************************************************/

/**
 * Data recorded by the SoftwareEventsProfiler
 */
struct SoftwareEventsSnapshot{
    uint64_t m_page_faults = 0; // is this perf::PERF_COUNT_SW_PAGE_FAULTS_MIN + perf::PERF_COUNT_SW_PAGE_FAULTS_MAJ ?
    uint64_t m_page_faults_min = 0; // requests handled by the O.S. page cache
    uint64_t m_page_faults_maj = 0; // Disk I/Os, the page is not in the O.S. cache
    uint64_t m_context_switches = 0; // because... why not?
    uint64_t m_cpu_migrations = 0; // not sure whether this is per thread or per process

    void operator+=(SoftwareEventsSnapshot snapshot);
    Database::BaseRecord data_record() const;
};

struct SoftwareEventsProfiler : public details::GenericProfiler {
    SoftwareEventsSnapshot m_current_snapshot;

public:
    /**
     * Initialise the profiler
     */
    SoftwareEventsProfiler();

    /**
     * Destructor
     */
    ~SoftwareEventsProfiler();


    /**
     * Start recording
     */
    void start();

    /**
     * Retrieve the data associated to this profiler
     */
    SoftwareEventsSnapshot snapshot();

    /**
     * Stop the recording
     */
    SoftwareEventsSnapshot stop();

    /**
     * Retrieve a data record ready to be stored in the database
     */
    Database::BaseRecord data_record();
};

std::ostream& operator<<(std::ostream& out, const SoftwareEventsSnapshot& snapshot);

} // namespace common


#endif //COMMON_PROFILER_HPP
