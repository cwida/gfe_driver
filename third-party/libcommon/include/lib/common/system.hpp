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

#ifndef COMMON_SYSTEM_HPP
#define COMMON_SYSTEM_HPP

#include <string>

namespace common {

/**
 * Retrieve the current hostname
 */
std::string hostname();

/**
 * Try to retrieve the last git commit for the current program (not libcommon). It returns an empty string in case of failure.
 * The idea is that a build is performed inside some directory from the source, e.g. <src>/build. Then the git repository
 * should be visible as part of <src>.
 */
std::string git_last_commit(); // impl in system_introspection.cpp

/**
 * Retrieve the amount of memory (RAM) available in the system, in bytes
 */
uint64_t get_total_ram();

/**
 * The information parsed from /proc/self/statm, apparently as multiples of 4Kb pages
 */
struct Statm {
    uint64_t m_vmsize; // total virtual memory, in 4kb pages
    uint64_t m_rss; // resident set size, in 4kb pages
    uint64_t m_shared; // number of resident shared pages (i.e., backed by a file)
    uint64_t m_text; // text segment, in 4kb pages
    uint64_t m_lib; // always 0
    uint64_t m_data; // data + stack, in 4kb pages
    uint64_t m_dt; // always 0
};

/**
 * Retrieve the content of /proc/self/statm
 */
Statm statm();

/**
 * Retrieve the total memory footprint of this process (excluding the text segment), in bytes
 */
uint64_t get_memory_footprint();

namespace filesystem {

/**
 * Retrieve the absolute path to the program executable
 */
std::string path_executable(); // impl in filesystem.cpp

/**
 * Retrieve the absolute path to the directory of the executable
 */
std::string directory_executable(); // impl in filesystem.cpp

} // namespace filesystem

/**
 * Concurrency related settings
 */
namespace concurrency {

/**
 * Get the Linux thread id, the value shown in the debugger
 */
int64_t get_thread_id();

/**
 * Get the Linux process id
 */
int64_t get_process_id();

/**
 * Get the number of threads alive in the current process
 */
int get_num_threads();

/**
 * Whether NUMA settings are available
 */
bool has_numa();

/**
 * Get the processor ID where the current thread is running
 */
int get_current_cpu();

/**
 * Get the Numa node associated to the CPU where the current thread is running
 */
int get_current_numa_node();

/**
 * Get the Numa node for the given CPU
 */
int get_numa_id(int cpu_id);

/**
 * Get the highest numa node in the system. Wrapper to `numa_max_node()'
 */
int get_numa_max_node();

/**
 * Pin the current thread to the current cpu and numa node. Optionally also pin the memory allocations from other NUMA nodes.
 */
void pin_thread_to_cpu(bool pin_numa_node = true);

/**
 * Pin the current thread to the given cpu and numa node. Disable memory allocations from the other NUMA nodes.
 * Note: as new created threads will inherit the same cpu mask, it is important to invoke this call to pin
 * an execution/sequential worker, rather than the main thread.
 */
void pin_thread_to_cpu(int cpu_id, bool pin_numa_node = true);

/**
 * Pin the current thread to the CPUs running at the given NUMA node
 */
void pin_thread_to_numa_node(int numa_node);

/**
 * Reset the pinning of the current thread
 * @param numa: if true, also unset the pinning to the NUMA node
 */
void unpin_thread(bool numa = true);

/**
 * Set the name of the current thread, it must a string up to 16 chars. The name is useful for debugging as
 * it is shown in the debugger's thread list.
 * @param name the new name to give to the current thread
 */
void set_thread_name(const std::string& name);

/**
 * Set the name of the given thread. The name must a string up to 16 chars. This is useful for debugging as
 * the thread's name is shown in the debugger's thread list.
 * @param thread_id the id the for the thread (pthread_t) which we want to set the name
 * @param name the new name to give to the given thread
 */
void set_thread_name(uint64_t thread_id, const std::string& name);

} // concurrency


/**
 * The compiler family being used
 */
enum class CompilerFamily { UNKNOWN, GCC, CLANG, INTEL };

/**
 * Detect on which compiler the source file has been compiled
 */
class CompilerInfo {
    CompilerFamily m_family = CompilerFamily::UNKNOWN; // the compiler family
    std::string m_name; // the name of the compiler
    int m_major = 0; // the major version
    int m_minor = 0; // the minor version
    int m_patch = 0; // the patch level

public:
    /**
     * Detect the compiler being used to build this file
     */
    CompilerInfo();

    /**
     * Get the compiler family
     */
    CompilerFamily family() const;

    /**
     * Get the name of the compiler
     */
    std::string name() const;

    /**
     * Get the major version of the compiler
     */
    int major() const;

    /**
     * Get the minor version of the compiler
     */
    int minor() const;

    /**
     * Get the patch level of the compiler
     */
    int patch() const;

    /**
     * Get a string description of the version of the compiler
     */
    std::string version() const;

    /**
     * Get a string representation of the compiler detected
     */
    std::string to_string() const;
};

// Overload of the out operator
std::ostream& operator<<(std::ostream& out, const CompilerInfo& ci);


/******************************************************************************
 *                                                                            *
 * Implementation details                                                     *
 *                                                                            *
 *****************************************************************************/
inline
CompilerInfo::CompilerInfo(){
#if defined(__clang__)
    m_family = CompilerFamily::CLANG;
    m_name = "Clang";
    m_major = __clang_major__;
    m_minor = __clang_minor__;
    m_patch = __clang_patchlevel__;
#elif defined(__INTEL_COMPILER)
    m_family = CompilerFamily::INTEL;
    m_name = "ICC";
    m_major = __INTEL_COMPILER / 100;
    m_minor = __INTEL_COMPILER % 100;
    m_patch = __INTEL_COMPILER_UPDATE;
#elif defined(__GNUC__)
    m_family = CompilerFamily::GCC;
    m_name = "GCC";
    m_major = __GNUC__;
    m_minor = __GNUC_MINOR__;
    m_patch = __GNUC_PATCHLEVEL__;
#endif
}


} // common


#endif //COMMON_SYSTEM_HPP
