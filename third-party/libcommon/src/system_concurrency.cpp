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

#include "system.hpp"

#include <cassert>
#include <cerrno>
#include <cstdio> // popen
#include <cstring> // strerror
#include <memory>
#include <pthread.h>
#include <sched.h> // sched_getcpu()
#include <sys/syscall.h> // SYS_gettid
#include <sys/sysinfo.h> // get_nprocs(): return the number of available processors
#include <unistd.h> // syscall

#include "error.hpp"

// Support for libnuma
#if __has_include(<numa.h>)
#include <numa.h>
#if !defined(HAVE_LIBNUMA)
#define HAVE_LIBNUMA
#endif
#endif

using namespace std;

namespace common::concurrency {

int64_t get_thread_id(){
    auto tid = (int64_t) syscall(SYS_gettid);
    assert(tid > 0);
    return tid;
}

int64_t get_process_id(){
    auto pid = (int64_t) syscall(SYS_getpid);
    assert(pid > 0);
    return pid;
}


int get_num_threads(){
    stringstream ss;
    ss << "cat /proc/" << get_process_id() << "/status | awk '/Threads/ {print $2}'";
    string stmt = ss.str();
    FILE* file = popen(stmt.c_str(), "r");
    if(file == nullptr) ERROR("Cannot retrieve the number of threads: " << strerror(errno) << " (errno: " << errno << ")");
    int num_threads = 0;
    if (fscanf(file, "%d", &num_threads) == EOF){
        num_threads = -1;
    }
    pclose(file);

    if(num_threads <= 0){ ERROR("Cannot retrieve the number of threads"); }

    return num_threads;
}

bool has_numa(){
#if defined(HAVE_LIBNUMA)
    return numa_available() != -1;
#else
    return false;
#endif
}

int get_current_cpu(){
    return sched_getcpu();
}

int get_current_numa_node(){
    return get_numa_id(get_current_cpu());
}

int get_numa_id(int cpu_id){
    int numa_node = -1;
#if defined(HAVE_LIBNUMA)
    if( has_numa() )
        numa_node = numa_node_of_cpu(cpu_id);
#endif
    return numa_node;
}

int get_numa_max_node(){
    int max_node = -1;
#if defined(HAVE_LIBNUMA)
    if( has_numa() )
        max_node = numa_max_node();
#endif
    return max_node;
}

void pin_thread_to_cpu(bool pin_numa_node){
    pin_thread_to_cpu(get_current_cpu(), pin_numa_node);
}

void pin_thread_to_cpu(int target_cpu, bool pin_numa_node){
    auto current_thread = pthread_self();
    cpu_set_t cpu_set;
    CPU_ZERO(&cpu_set);
    CPU_SET(target_cpu, &cpu_set);
    auto rc = pthread_setaffinity_np(current_thread, sizeof(cpu_set), &cpu_set);
    if (rc != 0) { ERROR("[pin_thread_to_cpu] pthread_setaffinity_np, rc: " << rc); }

#if defined(HAVE_LIBNUMA)
    if(pin_numa_node){
        const bool use_libnuma { has_numa() };

        if (use_libnuma) {
            auto current_node = use_libnuma ? numa_node_of_cpu(target_cpu) : -1;
//          numa_set_bind_policy(true); // not required: numa_set_membind is already strict (MPOL_BIND)
            auto nodemask_deleter = [](struct bitmask* ptr) {
                numa_free_nodemask(ptr);
            };
            unique_ptr<struct bitmask, decltype(nodemask_deleter)> nodemask_ptr {
                    numa_allocate_nodemask(), nodemask_deleter };
            auto nodemask = nodemask_ptr.get();
            numa_bitmask_setbit(nodemask, current_node);
            numa_set_membind(nodemask);
        }
    }
#endif
}

static void throw_error_numa_not_available(){
    ERROR("[pin_thread_to_numa_node] NUMA is not available in this system");
}

void pin_thread_to_numa_node(int numa_node){
    if(!has_numa()) throw_error_numa_not_available();
#if defined(HAVE_LIBNUMA)
    int rc = numa_run_on_node(numa_node);
    if(rc != 0){
       ERROR("[pin_thread_to_numa_node] Cannot pin the given node: " << numa_node << ", rc: " << rc << ", error: " << strerror(errno) << " (" << errno << "), ");
    }
#endif
}

void unpin_thread(bool unpin_numa) {
    auto current_thread = pthread_self();

    cpu_set_t cpu_set;
    CPU_ZERO(&cpu_set);
    for(size_t i = 0, n = get_nprocs(); i < n; i++){ CPU_SET(i, &cpu_set); }
    auto rc = pthread_setaffinity_np(current_thread, sizeof(cpu_set), &cpu_set);
    if (rc != 0) { ERROR("[unpin_thread] pthread_setaffinity_np, rc: " << rc); }

#if defined(HAVE_LIBNUMA)
    if( unpin_numa && has_numa() ){
        numa_set_localalloc(); // MPOL_DEFAULT
    }
#endif
}

void set_thread_name(const std::string& name){
    set_thread_name(pthread_self(), name);
}

void set_thread_name(uint64_t thread_id, const std::string& name) {
    std::string truncated_name = name.substr(0, 15);
    int rc = pthread_setname_np(thread_id, truncated_name.c_str());
    if (rc != 0) {ERROR("[set_thread_name] error: " << strerror(errno) << " (" << errno << ")"); }
}

} // common::concurrency
