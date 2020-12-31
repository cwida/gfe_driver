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

#include "memory_usage.hpp"

#include <cerrno>
#include <cstdlib>
#include <cstring> // strerror
#include <dlfcn.h> // dlsym
#include <elf.h>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>

#include "common/error.hpp"
#include "common/filesystem.hpp"
#include "common/system.hpp"

using namespace std;

// External variable, defined by POSIX, with the array of the env vars of the program
extern char **environ;

// Default name of the library to preload. It can be overloaded with the env. var. GFE_MEMORY_PROFILER
static const char* DEFAULT_LIBRARY_NAME_MEMORY_PROFILER = "gfe_memory_profiler.so";

// Pointer to the actual routine to compute the memory footprint
static int64_t (*fn_compute_memory_footprint) () = nullptr;

namespace gfe::utility {

static string path_memory_profiler(){
    // the name of the library to load
    const char* library_name = getenv("GFE_MEMORY_PROFILER");
    if(library_name == nullptr){ library_name = DEFAULT_LIBRARY_NAME_MEMORY_PROFILER; }

    // absolute path of the library
    string path;

    // search the path as given
    if(common::filesystem::file_exists(library_name)){
        return common::filesystem::absolute_path(library_name);
    }

    // search among the paths in LD_LIBRARY_PATH
    const char* ld_library_path = getenv("LD_LIBRARY_PATH");
    if(ld_library_path != nullptr){
        string directory;
        istringstream iss { ld_library_path };
        while(getline(iss, directory, ':')){
            path = directory + "/" + library_name;
            if(common::filesystem::file_exists(path)){
                return common::filesystem::absolute_path(path);
            }
        }
    }

    // search inside the executable directory
    path = common::filesystem::directory_executable() + "/" + library_name;
    if(common::filesystem::file_exists(path)){
        return common::filesystem::absolute_path(path);
    }

    // search in /lib
    path = string{"/lib/"} + library_name;
    if(common::filesystem::file_exists(path)){
        return common::filesystem::absolute_path(path);
    }

    // search in /usr/lib
    path = string{"/usr/lib/"} + library_name;
    if(common::filesystem::file_exists(path)){
        return common::filesystem::absolute_path(path);
    }

    ERROR("Unable to locate the library for the memory profiler: " << library_name);
}

void MemoryUsage::initialise(int argc, char* argv[]){
    if(is_initialised()) return; // already initialised

    // inform the loader of the library to override #malloc/#free
    string library_path = path_memory_profiler();
    setenv("LD_PRELOAD", library_path.c_str(), /* overwrite ? */ true);

    // reload the program
    string path_executable = common::filesystem::path_executable();
    execve(path_executable.c_str(), argv, environ);

    // we should never reach this point
    ERROR("Cannot reload the program `" << path_executable << ": " << strerror(errno)); // errno set by execve
}

bool MemoryUsage::is_initialised() {
    if(fn_compute_memory_footprint == nullptr){
        fn_compute_memory_footprint = reinterpret_cast<decltype(fn_compute_memory_footprint)>(dlsym(RTLD_DEFAULT, "gfe_compute_memory_footprint"));
    }

    return fn_compute_memory_footprint != nullptr;
}

int64_t MemoryUsage::memory_footprint(){
    if(fn_compute_memory_footprint != nullptr){ // trampoline to the actual implementation
        return fn_compute_memory_footprint();
    } else {
        return 0;
    }
}

uint64_t MemoryUsage::get_allocated_space(const void* pointer){
    if(pointer == nullptr) return 0;
    return reinterpret_cast<const uint64_t*>(pointer)[-1] & /* glibc flags */ ~7ull;
}

} // namespace


