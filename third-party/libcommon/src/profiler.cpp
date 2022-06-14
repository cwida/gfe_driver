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

#include "profiler.hpp"

#include <cassert>
#include <cstring>
#include <iostream>
#include <mutex>
#include <papi.h>

#undef CURRENT_ERROR_TYPE
#define CURRENT_ERROR_TYPE ::common::ProfilerError

using namespace std;

/*****************************************************************************
 *                                                                           *
 *   BaseProfiler                                                            *
 *                                                                           *
 *****************************************************************************/
namespace common { namespace details {

static std::mutex g_mutex;
bool BaseProfiler::library_initialised = false;

void BaseProfiler::initialise_library(){
    scoped_lock<mutex> lock(g_mutex);
    if(library_initialised) return; // already initialised
    int rc = PAPI_library_init(PAPI_VER_CURRENT);
    if (rc != PAPI_VER_CURRENT){ ERROR("Library PAPI version mismatch"); }
    library_initialised = true;
}

BaseProfiler::BaseProfiler() { initialise_library(); }

int BaseProfiler::get_event_code(const char* event_name){
    scoped_lock<mutex> lock(g_mutex);
    if(!library_initialised) ERROR("Library PAPI not initialised");

    // PAPI_event_name_to_code doesn't accept a const char* argument
    char buffer[PAPI_MAX_STR_LEN];
    strncpy(buffer, event_name, PAPI_MAX_STR_LEN);
    buffer[PAPI_MAX_STR_LEN -1] = '\0';
    int event_code { 0 };
    int rc = PAPI_event_name_to_code(buffer, &event_code);
    if(rc != PAPI_OK){ return -1; }

    PAPI_event_info_t event_info;
    rc = PAPI_get_event_info(event_code, &event_info);
    if(rc != PAPI_OK){ return -1; }

    if(IS_PRESET(event_info.event_code)){
        // check whether there exist some derivations to infer this event
        if(event_info.count > 0){
            return event_code;
        } else {
            return -1;
        }
    } else { // native event
        return event_code;
    }
}

}} // common::details

/*****************************************************************************
 *                                                                           *
 *   GenericProfiler                                                         *
 *                                                                           *
 *****************************************************************************/
namespace common { namespace details {

GenericProfiler::GenericProfiler(){ }

GenericProfiler::~GenericProfiler(){
    if(m_event_set != PAPI_NULL)
        unregister_events();
}

void GenericProfiler::add_events(const char* errorstring, const char* event_name){
    add_events(errorstring, &event_name, 1);
}

void GenericProfiler::add_events(const char* errorstring, const char* alternative_events[], size_t num_alternative_events){
    if (m_events_sz >= m_events_capacity) ERROR("No space left to add the events: " << m_events_sz);
    int rc = -1;
    size_t i = 0;
    while(rc == -1 && i < num_alternative_events){
        rc = get_event_code(alternative_events[i]);

        // add the event code && stop the execution
        if(rc != -1){
            m_events[m_events_sz] = rc;
            m_events_sz++;
        }

        i++;
    }

    if(rc == -1) {
        cerr << "[" __FILE__ << ":" << __LINE__ << "] " << errorstring << endl;
        ERROR(errorstring);
    }
}

void GenericProfiler::register_events(){
    scoped_lock<mutex> lock(g_mutex); // global access to PAPI
    int rc {0};

    m_event_set = PAPI_NULL;
    rc = PAPI_create_eventset(&m_event_set);
    if(rc != PAPI_OK) {
        cerr << "[" __FILE__ << ":" << __LINE__ << "] PAPI_create_eventset: " << PAPI_strerror(rc) << " (" << rc << ")" << endl;
        ERROR("Cannot create the event set (opaque object identifier for the PAPI library)");
    }

    rc = PAPI_add_events(m_event_set, m_events, m_events_sz);
    if(rc != PAPI_OK) {
        cerr << "[" __FILE__ << ":" << __LINE__ << "] PAPI_add_events: " << PAPI_strerror(rc) << " (" << rc << ")" << endl;
        ERROR("Cannot trace the interested set of events in this architecture");
    }
}

void GenericProfiler::unregister_events(){
    scoped_lock<mutex> lock(g_mutex);  // global access to PAPI
    int event_state = 0;
    int rc = PAPI_state(m_event_set, &event_state);
    if(rc == PAPI_OK && !(event_state & PAPI_STOPPED)){
        PAPI_stop(m_event_set, nullptr); // ignore rc
    }

    rc = PAPI_remove_events(m_event_set, m_events, m_events_sz);
    if(rc != PAPI_OK)
        cerr << "[" __FILE__ << ":" << __LINE__ << "] PAPI_remove_events: " << PAPI_strerror(rc) << " (" << rc << ")" << endl;

    rc = PAPI_destroy_eventset(&m_event_set);
    if(rc != PAPI_OK)
        cerr << "[" __FILE__ << ":" << __LINE__ << "] PAPI_destroy_eventset: " << PAPI_strerror(rc) << " (" << rc << ")" << endl;
}

void GenericProfiler::start(){
    int rc = PAPI_start(m_event_set);
    if(rc != PAPI_OK){
        ERROR("[GenericProfiler::start] Cannot start the event set: " << PAPI_strerror(rc) << "(rc: " << rc << ")");
    }
}

void GenericProfiler::stop(long long* resultset){
    int rc = PAPI_stop(m_event_set, (long long*) resultset);
    if(rc != PAPI_OK){
        ERROR("[GenericProfiler::stop] Cannot stop the event set: " << PAPI_strerror(rc) << "(rc: " << rc << ")");
    }
}

void GenericProfiler::snapshot(long long* resultset){
    int rc = PAPI_accum(m_event_set, resultset);
    if(rc != PAPI_OK){
        ERROR("[GenericProfiler::snapshot] Cannot obtain a snapshot from the event set: " << PAPI_strerror(rc) << "(rc: " << rc << ")");
    }
}

}} // common::details

/*****************************************************************************
 *                                                                           *
 *   CachesProfiler                                                          *
 *                                                                           *
 *****************************************************************************/
namespace common {

CachesProfiler::CachesProfiler() {
    add_events("Cannot infer cache-1 faults", "PAPI_L1_DCM");
    // on my damn AMD box L3 events are uncore :/
    const char* LLC_events[] = {"PAPI_L3_DCM", "PAPI_L3_TCM", "LLC-LOAD-MISSES"};
    add_events("Cannot infer cache-3 faults", LLC_events, 3);
    const char* TLB_events[] = {"PAPI_TLB_DM", "PAPI_TLB_TM"};
    add_events("Cannot infer TLB misses", TLB_events, 2);
    register_events();
}

CachesProfiler::~CachesProfiler() { }

void CachesProfiler::start(){
    GenericProfiler::start();
}

CachesSnapshot CachesProfiler::snapshot(){
    static_assert((sizeof(long long) * 3) == sizeof(m_current_snapshot), "Size mismatch, need to pass an array of types `long long'");
    GenericProfiler::snapshot((long long*) &m_current_snapshot);
    return m_current_snapshot;
}

CachesSnapshot CachesProfiler::stop(){
    CachesSnapshot m_result;

    static_assert(sizeof(long long) *3 == sizeof(m_result), "Size mismatch, need to pass an array of types `long long'");
    GenericProfiler::stop((long long*) &m_result);
    m_result += m_current_snapshot;

    m_current_snapshot = {0,0,0};

    return m_result;
}

Database::BaseRecord CachesProfiler::data_record(){
    return snapshot().data_record();
}

void CachesSnapshot::operator+=(CachesSnapshot snapshot){
    m_cache_l1_misses += snapshot.m_cache_l1_misses;
    m_cache_llc_misses += snapshot.m_cache_llc_misses;
    m_cache_tlb_misses += snapshot.m_cache_tlb_misses;
}

Database::BaseRecord CachesSnapshot::data_record() const {
    Database::BaseRecord record;
    record.add("cache_l1_misses", m_cache_l1_misses);
    record.add("cache_llc_misses", m_cache_llc_misses);
    record.add("cache_tlb_misses", m_cache_tlb_misses);
    return record;
}

std::ostream& operator<<(std::ostream& out, const CachesSnapshot& snapshot){
    out << "L1 faults: " << snapshot.m_cache_l1_misses << ", " <<
        "LLC faults: " << snapshot.m_cache_llc_misses << ", " <<
        "TLB faults: " << snapshot.m_cache_tlb_misses;
    return out;
}


CachesSnapshot operator+(const CachesSnapshot& s1, const CachesSnapshot& s2){
    CachesSnapshot result = s1;
    result.m_cache_l1_misses += s2.m_cache_l1_misses;
    result.m_cache_llc_misses += s2.m_cache_llc_misses;
    result.m_cache_tlb_misses += s2.m_cache_tlb_misses;
    return result;
}

/*****************************************************************************
 *                                                                           *
 *   BranchMispredictionsProfiler                                            *
 *                                                                           *
 *****************************************************************************/

BranchMispredictionsProfiler::BranchMispredictionsProfiler() {
    add_events("Cannot infer conditional branches", "PAPI_BR_CN");
    add_events("Cannot infer branch mispredictions", "PAPI_BR_MSP");
    // on my damn AMD box L3 events are uncore :/
    add_events("Cannot infer cache-1 faults", "PAPI_L1_DCM");
    const char* LLC_events[] = {"PAPI_L3_DCM", "PAPI_L3_TCM", "LLC-LOAD-MISSES"};
    add_events("Cannot infer cache-3 faults", LLC_events, 3);
    register_events();
}

BranchMispredictionsProfiler::~BranchMispredictionsProfiler() { }

void BranchMispredictionsProfiler::start() {
    GenericProfiler::start();
}

BranchMispredictionsSnapshot BranchMispredictionsProfiler::snapshot() {
    static_assert((sizeof(long long) * 4) == sizeof(m_current_snapshot), "Size mismatch, need to pass an array of types `long long'");
    GenericProfiler::snapshot((long long*) &m_current_snapshot);
    return m_current_snapshot;
}

BranchMispredictionsSnapshot BranchMispredictionsProfiler::stop() {
    BranchMispredictionsSnapshot m_result;

    static_assert(sizeof(long long) *4 == sizeof(m_result), "Size mismatch, need to pass an array of types `long long'");
    GenericProfiler::stop((long long*) &m_result);
    m_result += m_current_snapshot;

    m_current_snapshot = {0,0,0};

    return m_result;
}

Database::BaseRecord BranchMispredictionsProfiler::data_record(){
    return snapshot().data_record();
}


void BranchMispredictionsSnapshot::operator+=(BranchMispredictionsSnapshot snapshot) {
    m_conditional_branches += snapshot.m_conditional_branches;
    m_branch_mispredictions += snapshot.m_branch_mispredictions;
    m_cache_l1_misses += snapshot.m_cache_l1_misses;
    m_cache_llc_misses += snapshot.m_cache_llc_misses;
}

Database::BaseRecord BranchMispredictionsSnapshot::data_record() const {
    Database::BaseRecord record;
    record.add("conditional_branches", m_conditional_branches);
    record.add("branch_mispredictions", m_branch_mispredictions);
    record.add("cache_l1_misses", m_cache_l1_misses);
    record.add("cache_llc_misses", m_cache_llc_misses);
    return record;
}

std::ostream& operator<<(std::ostream& out, const BranchMispredictionsSnapshot& snapshot){
    out << "Conditional branches: " << snapshot.m_conditional_branches << ", " <<
        "Branch mispredictions: " << snapshot.m_branch_mispredictions << ", " <<
        "L1 cache faults: " << snapshot.m_cache_l1_misses << ", " <<
        "LLC cache faults: " << snapshot.m_cache_llc_misses;
    return out;
}

/*****************************************************************************
 *                                                                           *
 *   SoftwareEventsProfiler                                                  *
 *                                                                           *
 *****************************************************************************/
SoftwareEventsProfiler::SoftwareEventsProfiler() {
    add_events("Cannot install PERF_COUNT_SW_PAGE_FAULTS", "perf::PERF_COUNT_SW_PAGE_FAULTS");
    add_events("Cannot install PERF_COUNT_SW_PAGE_FAULTS_MIN", "perf::PERF_COUNT_SW_PAGE_FAULTS_MIN");
    add_events("Cannot install PERF_COUNT_SW_PAGE_FAULTS_MAJ", "perf::PERF_COUNT_SW_PAGE_FAULTS_MAJ");
    add_events("Cannot install PERF_COUNT_SW_CONTEXT_SWITCHES", "perf::PERF_COUNT_SW_CONTEXT_SWITCHES");
    add_events("Cannot install PERF_COUNT_SW_CPU_MIGRATIONS", "perf::PERF_COUNT_SW_CPU_MIGRATIONS");
    register_events();
}

SoftwareEventsProfiler::~SoftwareEventsProfiler() { }

void SoftwareEventsProfiler::start() {
    GenericProfiler::start();
}

SoftwareEventsSnapshot SoftwareEventsProfiler::snapshot() {
    static_assert((sizeof(long long) * 5) == sizeof(m_current_snapshot), "Size mismatch, need to pass an array of types `long long'");
    GenericProfiler::snapshot((long long*) &m_current_snapshot);
    return m_current_snapshot;
}

SoftwareEventsSnapshot SoftwareEventsProfiler::stop() {
    SoftwareEventsSnapshot m_result;

    static_assert(sizeof(long long) *5 == sizeof(m_result), "Size mismatch, need to pass an array of types `long long'");
    GenericProfiler::stop((long long*) &m_result);
    m_result += m_current_snapshot;

    m_current_snapshot = {0,0,0,0,0};

    return m_result;
}

Database::BaseRecord SoftwareEventsProfiler::data_record(){
    return snapshot().data_record();
}

void SoftwareEventsSnapshot::operator+=(SoftwareEventsSnapshot snapshot) {
    m_page_faults += snapshot.m_page_faults;
    m_page_faults_min += snapshot.m_page_faults_min;
    m_page_faults_maj += snapshot.m_page_faults_maj;
    m_context_switches += snapshot.m_context_switches;
    m_cpu_migrations += snapshot.m_cpu_migrations;
}

Database::BaseRecord SoftwareEventsSnapshot::data_record() const {
    Database::BaseRecord record;
    record.add("page_faults", m_page_faults);
    record.add("page_faults_minor", m_page_faults_min);
    record.add("page_faults_major", m_page_faults_maj);
    record.add("context_switches", m_context_switches);
    record.add("cpu_migrations", m_cpu_migrations);
    return record;
}

std::ostream& operator<<(std::ostream& out, const SoftwareEventsSnapshot& snapshot){
    out << "Page faults: " << snapshot.m_page_faults << " "
           "(minor: " << snapshot.m_page_faults_min << ", major: " << snapshot.m_page_faults_maj << "), " <<
           "Context switches: " << snapshot.m_context_switches << ", "
           "CPU migrations: " << snapshot.m_cpu_migrations;
    return out;
}

} // namespace common

