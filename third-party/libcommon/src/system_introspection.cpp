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
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <sys/sysinfo.h> // linux specific
#include <unistd.h> // chdir, gethostname

#include "error.hpp"
#include "filesystem.hpp"

using namespace std;

namespace common {

string hostname(){
    constexpr int len = 512;
    char buffer[len];
    auto rc = gethostname(buffer, len);
    if(rc != 0){
        ERROR("[hostname] Cannot retrieve the hostname: " << strerror(errno)  << " (" << errno << ")");
    }
    string hostname{buffer};

    // Remove the suffix `.scilens.private' from the machines in the Scilens cluster
    const string scilens_suffix = ".scilens.private";
    if(hostname.length() >= scilens_suffix.length()){
        auto scilens_match = hostname.rfind(scilens_suffix);
        if(scilens_match != string::npos && scilens_match == hostname.length() - scilens_suffix.length()){
            hostname.replace(scilens_match, scilens_suffix.length(), nullptr, 0);
        }
    }

    return hostname;
}

// Helper, find the source directory from CMake
static string git_get_srcdir_from_cmake(){
    if(!filesystem::exists("CMakeCache.txt")) return "";
    string result;

    fstream f("CMakeCache.txt", ios::in);
    regex pattern("^CMAKE_HOME_DIRECTORY(:.+)?=(.+?)\\s*$");

    while(!f.eof()){
        string line;
        std::getline(f, line);
        smatch matches;
        if( regex_match(line, matches, pattern) ){
            assert(matches.size() == 3);
            result = matches[2];
            break;
        }
    }
    f.close();

    return result;
}


// Helper, find the source directory from the Makefile
static string git_get_srcdir_from_makefile(){
//    if(!fs::exists("Makefile")) return "";
    if(!filesystem::exists("Makefile")) return "";

    constexpr size_t buffer_sz = 512;
    char buffer[buffer_sz];

    FILE* fp = popen("make -n -p | awk '/^srcdir\\s*:?=/ {print $NF}'", "r");
    if(fp == nullptr){
        cerr << "[get_srcdir_from_makefile] WARNING: Cannot execute make: " << strerror(errno) << endl;
        return "";
    }
    char* result = fgets(buffer, buffer_sz, fp);
    if(result != buffer){
        cerr << "[get_srcdir_from_makefile] WARNING: Cannot read the result from make: " << strerror(errno) << endl;
    } else {
        // (chomp) truncate the string at the first '\n'
        strtok(result, "\n");
    }

    fclose(fp);
    return result;
}

// Helper, try to execute `git' and grep for the last commit
static string git_read_last_commit(){
    FILE* fp = popen("git log -1 | awk 'NR == 1 && $1 == \"commit\" {print $2}'", "r");

    constexpr size_t buffer_sz = 512;
    char buffer[buffer_sz];

    if(fp == nullptr){
        cerr << "[git_read_last_commit] WARNING: Cannot retrieve the last git version: " << strerror(errno) << endl;
        return "";
    }

    char* result = fgets(buffer, buffer_sz, fp);
    if(result != buffer){
        cerr << "[git_read_last_commit] WARNING: Cannot read the result from git: " << strerror(errno) << endl;
    } else {
        // (chomp) truncate the string at the first '\n'
        strtok(result, "\n");
    }

    fclose(fp);
    return result;
}

string git_last_commit(){
    auto basedir = filesystem::directory_executable();
    filesystem::TemporaryWorkingDirectory tmpwd{basedir}; // restore the previous wd on exit

    // move to the source directory if possible
    auto srcdir = git_get_srcdir_from_cmake(); // try with cmake
    if(srcdir.empty()){ // try again with Makefile
        srcdir = git_get_srcdir_from_makefile();
    }
    if(!srcdir.empty()) {
        int rc = chdir(srcdir.c_str());
        if (rc != 0) {
            cerr << "[git_last_commit] ERROR: cannot change the current working directory to " << srcdir << ": "
                 << strerror(errno) << endl;
        }
    }

    // finally try to execute `git' and fetch the last commit
    return git_read_last_commit();
}

uint64_t get_total_ram(){
    struct sysinfo info;
    int rc = sysinfo(&info);
    if(rc != 0) ERROR("[get_total_ram] Invocation to sysinfo() failed: " << strerror(errno) << " (errno: " << errno << ")");
    return static_cast<uint64_t>(info.totalram) * info.mem_unit;
}

Statm statm() {
    fstream f("/proc/self/statm", ios_base::in);
    if(!f.good()){ ERROR("[statm] Cannot access /proc/self/statm"); }

    Statm s;
    f >> s.m_vmsize;
    f >> s.m_rss; // amount of physical memory used
    f >> s.m_shared;
    f >> s.m_text;
    f >> s.m_lib;
    f >> s.m_data; // virtual memory (incl. non allocated memory mapped files)
    f >> s.m_dt;

    f.close();
    return s;
}

uint64_t get_memory_footprint() {
    return statm().m_rss * /* 4 Kb */ (1<<12);
}

} // namespace common
