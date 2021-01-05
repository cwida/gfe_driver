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

#include "teseo_openmp.hpp"

#include <cassert>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <omp.h>
#include <sched.h>

#include "teseo_driver.hpp"

#include "teseo.hpp"

using namespace std;

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG
namespace gfe { extern mutex _log_mutex [[maybe_unused]]; }
#define COUT_DEBUG_CLASSNAME ""
#define COUT_DEBUG_FORCE(msg) { std::scoped_lock<std::mutex> lock{::gfe::_log_mutex}; std::cout << "[" << COUT_DEBUG_CLASSNAME << "::" << __FUNCTION__ << "] " << msg << std::endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif

namespace gfe::library::teseo_driver_internal {

/*****************************************************************************
 *                                                                           *
 *  ThreadAffinity                                                           *
 *                                                                           *
 *****************************************************************************/
#undef COUT_DEBUG_CLASSNAME
#define COUT_DEBUG_CLASSNAME "ThreadAffinity"

// 16 and 129 below are arbitrary values, enough to represent the current processor sizes
using processors_t = int [16][129]; // up to 16 sockets, the first value is the number of processors in socket
static constexpr int max_num_sockets = sizeof(processors_t{}) / sizeof(processors_t{}[0]);
static constexpr int max_num_processors = sizeof(processors_t{}[0]) / sizeof(processors_t{}[0][0]) -1;

static void parse_cpuinfo(int& num_sockets, processors_t& processors){
    constexpr int buffer_sz = 1024;
    char buffer[buffer_sz];

    // init
    for(int i = 0; i < max_num_sockets; i++){
        processors[i][0] = 0;
    }
    num_sockets = 0;

    FILE* file = fopen("/proc/cpuinfo", "r");
    if(file == nullptr){
        fprintf(stderr, "[parse_cpuinfo] ERROR: cannot open /proc/cpuinfo, errno: %d, error: %s\n", errno, strerror(errno));
        ::exit(EXIT_FAILURE);
    }

    bool read_processor_id = true; // parse either `processor : <value>` or `physical id: <value>`

    int processor_id = -1;
    int socket_id = -1;

    while(fgets(buffer, buffer_sz, file) != nullptr){
        uint64_t buffer_len = strlen(buffer);
        uint64_t key_len = 0;
        while(key_len < buffer_len && buffer[key_len] != ':'){ key_len++; }
        if(key_len == buffer_len) continue; // no key read

        if(read_processor_id){ // read the processor
            if(key_len < strlen("processor")) continue; // next line

            if(strncmp("processor", buffer, strlen("processor")) == 0){
                char* value = buffer + key_len +1;
                processor_id = atoi(value);

                read_processor_id = false; // now parse the socket id
            }
        } else { // read the physical id
            if(key_len < strlen("physical id")) continue; // next line

            if(strncmp("physical id", buffer, strlen("physical id")) == 0){
                char* value = buffer + key_len +1;
                socket_id = atoi(value);

                //printf("processor: %d, socket: %d\n", processor_id, socket_id);
                assert(socket_id >= 0 && processor_id >= 0 && "[parse_cpuinfo] Incorrect reading");

                if(socket_id >= max_num_processors){
                    fprintf(stderr, "[parse_cpuinfo] ERROR: invalid socket ID: %d, we can support up to 16 sockets\n", socket_id);
                    exit(EXIT_FAILURE);
                }

                if(socket_id >= num_sockets){ // assume the sockets are numbered 0, 1, 2, ...
                    num_sockets = socket_id +1;
                }

                if(processors[socket_id][0] >= max_num_processors){
                    fprintf(stderr, "[parse_cpuinfo] ERROR: already reached the max number of processors (%d) for the socket %d\n", max_num_processors, socket_id);
                    exit(EXIT_FAILURE);
                }


                processors[socket_id][1 + processors[socket_id][0]] = processor_id;
                processors[socket_id][0]++;

                read_processor_id = true; // the next iteration parse the processor id
                processor_id = socket_id = -1;
            }
        }
    }


    if(!feof(file)){ // we expect to have read the whole file
        fprintf(stderr, "[parse_cpuinfo] ERROR: cannot parse /proc/cpuinfo, errno: %d, error: %s\n", errno, strerror(errno));
        exit(EXIT_FAILURE);
    }

    fclose(file); file = nullptr;
}

namespace { // anon
class NodeAffinityMasks {
    NodeAffinityMasks(const NodeAffinityMasks&) = delete;
    NodeAffinityMasks& operator=(const NodeAffinityMasks&) = delete;

    int m_num_nodes;
    cpu_set_t m_cpu_sets[max_num_sockets];

public:
    NodeAffinityMasks(){
        processors_t processors;
        parse_cpuinfo(m_num_nodes, processors);

        // init
        for(int i = 0; i < max_num_sockets; i++){
            CPU_ZERO(&m_cpu_sets[i]);
        }

        for(int i = 0; i < m_num_nodes; i++){
            int num_physical_threads = processors[i][0];
            int* physical_threads = &processors[i][1];

            for(int j = 0; j < num_physical_threads; j++){
                //printf("node: %d, physical thread: %d\n", i, physical_threads[j]);
                CPU_SET(physical_threads[j], &m_cpu_sets[i]);
            }
        }
    }

    int num_nodes() const {
        return m_num_nodes;
    }

    const cpu_set_t* cpuset(int node_id) const {
        assert(node_id < num_nodes() && "Invalid node ID");
        return &m_cpu_sets[node_id];
    }
};
} // anon namespace

static NodeAffinityMasks g_node_affinity_masks;

ThreadAffinity::ThreadAffinity(bool enabled) : m_is_enabled(enabled) {
    if(m_is_enabled){
        int thread_id = omp_get_thread_num();

        if(thread_id == 0){
            m_copy_to_restore = (cpu_set_t*) malloc(sizeof(cpu_set_t));
            if(m_copy_to_restore == nullptr) {
                fprintf(stderr, "[ThreadAffinity::ctor] ERROR: Not enough space to allocate a cpu_set_t\n");
                exit(EXIT_FAILURE);
            }
            int rc = sched_getaffinity(/* current thread */ 0, sizeof(cpu_set_t), m_copy_to_restore);
            if(rc != 0){
                fprintf(stderr, "[ThreadAffinity::ctor] ERROR: call to sched_getaffinity failed, errno: %d, strerror: %s\n", errno, strerror(errno));
                exit(EXIT_FAILURE);
            }
        }

        int socket_id = thread_id % g_node_affinity_masks.num_nodes();
        int rc = sched_setaffinity(/* this thread */ 0, sizeof(cpu_set_t), g_node_affinity_masks.cpuset(socket_id));
        if(rc != 0){
            fprintf(stderr, "[ThreadAffinity::ctor] ERROR: call to sched_setaffinity failed, errno: %d, strerror: %s\n", errno, strerror(errno));
            exit(EXIT_FAILURE);
        }
    }
}

ThreadAffinity::~ThreadAffinity(){
    if(is_enabled() && m_copy_to_restore != nullptr){
        int rc = sched_setaffinity(/* this thread */ 0,  sizeof(cpu_set_t), m_copy_to_restore);
        if(rc != 0){
            fprintf(stderr, "[ThreadAffinity::dtor] ERROR: Cannot restore the cpu_set_t, errno: %d, strerror: %s\n", errno, strerror(errno));
            exit(EXIT_FAILURE);
        }
        free(m_copy_to_restore); m_copy_to_restore = nullptr;
    }
}

bool ThreadAffinity::is_enabled() const noexcept {
    return m_is_enabled;
}


/*****************************************************************************
 *                                                                           *
 *  RegisterThread                                                           *
 *                                                                           *
 *****************************************************************************/
#undef COUT_DEBUG_CLASSNAME
#define COUT_DEBUG_CLASSNAME "RegisterThread"

RegisterThread::RegisterThread(::teseo::Teseo* teseo) : m_encoded_addr(reinterpret_cast<uint64_t>(teseo)){
    assert(teseo != nullptr);
    assert(teseo == this->teseo() && "Cannot decode the address of teseo");
    assert(omp_get_thread_num() == 0 && "Expected to be invoked in the master thread");

    teseo->register_thread();
    assert(!is_enabled() && "Invalid address");
    m_encoded_addr |= 0x1ull;
    assert(is_enabled() && "That's the whole purpose of the encoding");
}

RegisterThread::RegisterThread(const RegisterThread& rt) : m_encoded_addr(reinterpret_cast<uint64_t>(rt.teseo())){
    if(omp_get_thread_num() > 0 ){
        teseo()->register_thread();
        m_encoded_addr |= 0x1ull;
    }
}

RegisterThread::~RegisterThread(){
    if(is_enabled()){
        teseo()->unregister_thread();
    }
}

::teseo::Teseo* RegisterThread::teseo() const {
    return reinterpret_cast<::teseo::Teseo*>(m_encoded_addr & ~0x1ull);
}

bool RegisterThread::is_enabled() const {
    return m_encoded_addr & 0x1ull;
}


/*****************************************************************************
 *                                                                           *
 *  OpenMP                                                                   *
 *                                                                           *
 *****************************************************************************/
#undef COUT_DEBUG_CLASSNAME
#define COUT_DEBUG_CLASSNAME "OpenMP"

OpenMP::OpenMP(TeseoDriver* driver) :
        m_thread_affinity(driver->has_thread_affinity()),
        m_thread_tracking(reinterpret_cast<teseo::Teseo*>(driver->handle_impl())),
        m_transaction(reinterpret_cast<teseo::Teseo*>(driver->handle_impl())->start_transaction(driver->has_read_only_transactions())),
        m_iterator(m_transaction.iterator()) {
    COUT_DEBUG("master ctor, omp thread: " << omp_get_thread_num() << ", thread_id: " << pthread_self());
}

OpenMP::OpenMP(const OpenMP& copy) :
        m_thread_affinity(copy.m_thread_affinity.is_enabled() && omp_get_thread_num() != 0), /* do not re-init the master thread */
        m_thread_tracking(copy.m_thread_tracking),
        m_transaction(copy.m_transaction),
        m_iterator(m_transaction.iterator()) {
    COUT_DEBUG("copy ctor, omp thread: " << omp_get_thread_num() << ", thread_id: " << pthread_self());
}

OpenMP::~OpenMP(){
    /* nop */
}

} // namespace
