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

#include <atomic>
#include <cassert>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib> // unsetenv
#include <cstring>
#include <dlfcn.h>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <sys/mman.h>
#include <sys/syscall.h>
#include <unistd.h>

using namespace std;

/*****************************************************************************
 *                                                                           *
 *  Init                                                                     *
 *                                                                           *
 *****************************************************************************/
static bool g_initialised = false; // whether the hooks have been initialised
static mutex g_mutex;
static atomic<int64_t> g_next_thread_id = 0;
static thread_local int g_thread_id = -1;
static thread_local bool g_recursion = false; // whether we are calling the function from glibc (recursive call)
static char g_popen_stmt[256] = {0}; // the statement to read the process mappings, set by #initialise_once

static constexpr uint64_t num_entries = 1ull << 16;
struct {
    int64_t m_size; // sum of allocated/free memory in the current thread
} g_thread_local_entries[num_entries];
static uint64_t g_memory_mappings[num_entries];
static uint64_t g_num_memory_mappings = 0;

// real functions
static void* (*glibc_malloc)(size_t sz) = nullptr;
static void* (*glibc_calloc)(size_t num, size_t sz) = nullptr;
static void* (*glibc_realloc) (void* ptr, size_t size) = nullptr;
static void (*glibc_free)(void* ptr) = nullptr;
static int (*glibc_posix_memalign)(void **memptr, size_t alignment, size_t size) = nullptr;
static void* (*glibc_mmap)(void *addr, size_t length, int prot, int flags, int fd, off_t offset);
static int (*glibc_munmap)(void *addr, size_t length);

// internally malloc uses the last 3 bits of the allocated size as control bits : PREV_INUSE, IS_MAPPED, NON_MAIN_ARENA
// mask them to get the actual allocated size (including overhead).
static constexpr uint64_t MASK_MALLOC_SIZE = ~ static_cast<uint64_t>((0x1 | 0x2 | 0x4));

// dlsym invokes calloc. We initially resolve these invocations through a small memory pool
void* gfe_dummy_calloc(size_t num, size_t size){
    uint64_t requested_bytes = num * size;

    constexpr uint64_t buffer_sz = 1024;
    static char buffer[buffer_sz] = {0};
    static uint64_t end = 0; // next available

    uint64_t start = end;

    end += requested_bytes;
    if(end >= buffer_sz){ // overflow
        fprintf(stderr, "[memory profiler error] depleted the space in the initial memory pool\n");
        exit(EXIT_FAILURE);
    }

    return buffer + start;
}

void initialise_once(){
    if(g_thread_id < 0)// ensure that a thread_id is set
        g_thread_id = g_next_thread_id++;

    // check-lock-check again
    if (g_initialised) return;
    unique_lock<std::mutex> xlock(g_mutex);
    if (g_initialised) return;

    // resolve the hooks
    glibc_calloc = gfe_dummy_calloc;
    glibc_calloc = reinterpret_cast<decltype(glibc_calloc)>(dlsym(RTLD_NEXT, "calloc"));
    if(glibc_calloc == nullptr){
        fprintf(stderr, "[memory profiler error] cannot resolve the function calloc\n");
        exit(EXIT_FAILURE);
    }
    glibc_malloc = reinterpret_cast<decltype(glibc_malloc)>(dlsym(RTLD_NEXT, "malloc"));
    if(glibc_malloc == nullptr){
        fprintf(stderr, "[memory profiler error] cannot resolve the function malloc\n");
        exit(EXIT_FAILURE);
    }
    glibc_realloc = reinterpret_cast<decltype(glibc_realloc)>(dlsym(RTLD_NEXT, "realloc"));
    if(glibc_malloc == nullptr){
        fprintf(stderr, "[memory profiler error] cannot resolve the function realloc\n");
        exit(EXIT_FAILURE);
    }
    glibc_free = reinterpret_cast<decltype(glibc_free)>(dlsym(RTLD_NEXT, "free"));
    if(glibc_malloc == nullptr){
        fprintf(stderr, "[memory profiler error] cannot resolve the function free\n");
        exit(EXIT_FAILURE);
    }
    glibc_posix_memalign = reinterpret_cast<decltype(glibc_posix_memalign)>(dlsym(RTLD_NEXT, "posix_memalign"));
    if(glibc_posix_memalign == nullptr){
        fprintf(stderr, "[memory profiler error] cannot resolve the function posix_memalign\n");
        exit(EXIT_FAILURE);
    }
    glibc_mmap = reinterpret_cast<decltype(glibc_mmap)>(dlsym(RTLD_NEXT, "mmap"));
    if(glibc_mmap == nullptr){
        fprintf(stderr, "[memory profiler error] cannot resolve the function mmap\n");
        exit(EXIT_FAILURE);
    }
    glibc_munmap = reinterpret_cast<decltype(glibc_munmap)>(dlsym(RTLD_NEXT, "munmap"));
    if(glibc_munmap == nullptr){
        fprintf(stderr, "[memory profiler error] cannot resolve the function munmap\n");
        exit(EXIT_FAILURE);
    }

    memset(g_thread_local_entries, '\0', num_entries * sizeof(g_thread_local_entries[0]));

    unsetenv("LD_PRELOAD"); // reset the env. var. used to load this library

    { // create the statement for popen
        stringstream ss;
        ss << "/usr/bin/pmap ";
        ss << (int64_t) syscall(SYS_getpid);
        ss << " -xq";
        string str = ss.str();
        memcpy(g_popen_stmt, str.c_str(), str.size());
    }

    g_initialised = true;
}


namespace {
struct RecursionGuard {
    RecursionGuard() {
        assert(g_recursion == false && "Recursion");
        g_recursion = true;
    }

    ~RecursionGuard(){
        g_recursion = false;
    }
};
} // anonymous namespace

/*****************************************************************************
 *                                                                           *
 *  Run-time                                                                 *
 *                                                                           *
 *****************************************************************************/

static void handle_malloc(void* pointer, uint64_t requested_size){
    if(pointer == nullptr) return; // nop

    int64_t allocated_size = reinterpret_cast<uint64_t*>(pointer)[-1] & MASK_MALLOC_SIZE;
    g_thread_local_entries[g_thread_id].m_size += allocated_size;

    // printf("[malloc] pointer: %p, requested_size: %lu, allocated size: %lld\n", pointer, requested_size, allocated_size);
}

static void handle_free(void* pointer){
    if(pointer == nullptr) return; // nop

    int64_t allocated_size = reinterpret_cast<uint64_t*>(pointer)[-1] & MASK_MALLOC_SIZE;
    g_thread_local_entries[g_thread_id].m_size -= allocated_size;

    // printf("[free] pointer: %p, allocated size: %ld\n", pointer, allocated_size);
}

static void handle_mmap(void* pointer){
    unique_lock<mutex> xlock(g_mutex);

    if(g_num_memory_mappings >= num_entries){
        fprintf(stderr, "[memory profiler error] Too many memory mappings! (ARGH!)\n");
        exit(EXIT_FAILURE);
    }

    uint64_t address = reinterpret_cast<uint64_t>(pointer);
    int64_t i = g_num_memory_mappings;
    while(i > 0 && g_memory_mappings[i -1] > address){
        g_memory_mappings[i] = g_memory_mappings[i -1];
        i--;
    }
    g_memory_mappings[i] = address;

    g_num_memory_mappings++;
}

static void handle_munmap(void* pointer){
    if(pointer == nullptr) return; // nop
    unique_lock<mutex> xlock(g_mutex);

    // well, we assume that the whole mapped region is unmapped in one go
    uint64_t address = reinterpret_cast<uint64_t>(pointer);

    int64_t i = 0, num_memory_mappings = g_num_memory_mappings;
    bool found = false;
    while(!found && i < num_memory_mappings){
        assert((i == 0 || g_memory_mappings[i] > g_memory_mappings[i-1]) && "Sorted order not respected");

        if(g_memory_mappings[i] == address){
            found = true;
        } else {
            i++;
        }
    }

    if(!found){
        fprintf(stderr, "[memory profiler error] memory mapping not found!");
        exit(EXIT_FAILURE);
    }

   while(i < num_memory_mappings){
       g_memory_mappings[i] = g_memory_mappings[i+1];
       i++;
   }

   g_num_memory_mappings--;
}

/*****************************************************************************
 *                                                                           *
 *  Hooks                                                                    *
 *                                                                           *
 *****************************************************************************/

#define PREAMBLE(action) \
    if(g_recursion){ return action; } /* avoid that the hooks being fired when the functions are invoked by glibc */ \
    RecursionGuard _recursion_guard; \
    initialise_once(); \


extern "C" {

void* malloc(size_t size) {
    PREAMBLE( glibc_malloc(size) );

    void* ret = glibc_malloc(size);
    handle_malloc(ret, size);
    return ret;
}

void* calloc(size_t num, size_t sz) {
    PREAMBLE( glibc_calloc(num, sz) );

    void* ret = glibc_calloc(num, sz);
    handle_malloc(ret, num * sz);
    return ret;
}

void* realloc(void* ptr, size_t size){
    PREAMBLE( glibc_realloc(ptr, size) );

    handle_free(ptr);
    void* ret = glibc_realloc(ptr, size);
    handle_malloc(ptr, size);

    return ret;
}

void free(void* ptr){
    PREAMBLE( glibc_free(ptr) );

    handle_free(ptr);
    glibc_free(ptr);
}

int posix_memalign(void **memptr, size_t alignment, size_t size){
    PREAMBLE( glibc_posix_memalign(memptr, alignment, size) );

    int rc = posix_memalign(memptr, alignment, size);
    if(rc == 0){ // success
        handle_malloc(*memptr, size);
    }

    return rc;
}

void* mmap(void *addr, size_t length, int prot, int flags, int fd, off_t offset) {
    PREAMBLE( glibc_mmap(addr, length, prot, flags, fd, offset) );

    void* ret = glibc_mmap(addr, length, prot, flags, fd, offset);
    if(ret != MAP_FAILED){ // record the mapping
        //printf("mmap, return address: %p (%" PRIu64 "), length: %zu, prot: %d, flags: %d, fd: %d, offset: %ld \n", ret, (uint64_t) ret, length, prot, flags, fd, offset);
        handle_mmap(ret);
    }

    return ret;
}

int munmap(void *addr, size_t length){
    PREAMBLE( glibc_munmap(addr, length) );


    int rc = glibc_munmap(addr, length);
    if(rc == 0){
        //printf("munmap, address: %p, length: %zu\n", addr, length);
        handle_munmap(addr);
    }

    return rc;
}

} // avoid mangling

/*****************************************************************************
 *                                                                           *
 *  API                                                                      *
 *                                                                           *
 *****************************************************************************/
extern "C" { // avoid mangling

int64_t gfe_compute_memory_footprint (){
    unique_lock<mutex> xlock(g_mutex);
    int64_t memfp = 0;

    // count the virtual memory allocated by glibc
    int64_t num_entries = g_next_thread_id;
    for(int64_t i = 0; i < num_entries; i++){
        memfp += g_thread_local_entries[i].m_size;
    }

    // count the amount of physical memory used by the memory mappings
    FILE* fp = popen(g_popen_stmt, "r");
    if(fp == nullptr){
        fprintf(stderr, "[memory profiler error] popen failed, stmt: %s\n", g_popen_stmt);
        exit(EXIT_FAILURE);
    }
    char buffer[1024];
    uint64_t g_next_mapping = 0;
    while(fgets(buffer, sizeof(buffer), fp) != nullptr && g_next_mapping < g_num_memory_mappings){
        //printf("%s", buffer); // it already terminates with a \n
        char* next = nullptr;
        uint64_t address = strtoull(buffer, &next, 16);
        if(address < g_memory_mappings[g_next_mapping]) continue; // skip
        while(g_next_mapping < g_num_memory_mappings && address > g_memory_mappings[g_next_mapping]){
            fprintf(stderr, "[memory profiler error] address not matched: %" PRIu64 "\n", g_memory_mappings[g_next_mapping]);
            g_next_mapping++;
        }
        if(g_next_mapping < g_num_memory_mappings && address == g_memory_mappings[g_next_mapping]){ // match
            strtoull(next, &next, 10); // skip the virtual memory
            uint64_t rss = strtoull(next, nullptr, 10) * /* KB */ 1024ull;
            //printf("address: %" PRIu64 ", rss: %" PRIu64 "\n", address, rss);

            g_next_mapping++;
            memfp += rss;
        }
    }
    pclose(fp);

    return memfp;
}

} // extern "C"
