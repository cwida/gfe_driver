/**
 * This file is part of libcommon.
 *
 * libcommon is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
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
#include "backtrace.hpp"

#include <cassert>
#include <cxxabi.h> // demangler
#include <cstdlib>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <stdexcept>
#include <thread>
#include "backtrace.h" // API libbacktrace

using namespace std;

namespace common {

/**
 * Global state
 */
static struct backtrace_state* g_backtrace_state = nullptr;
static mutex g_mutex;


/**
 * Callbacks
 */

static void on_error(void* instance, const char* msg, int errnum) {
    stringstream ss;
    ss << msg << "(" << errnum << ")";
    throw std::runtime_error(ss.str());
}

extern "C" {
int _on_backtrace_full(void* ptr_instance, uintptr_t pc, const char* filename, int lineno, const char* function) {
    assert(ptr_instance != nullptr);
    Backtrace* instance = reinterpret_cast<Backtrace*>(ptr_instance);
    if (instance->m_skip_next_frames) return 0; // avoid recording this frame

    string str_filename, function_mangled, function_demangled;
    if (filename != nullptr) {
        str_filename = filename;
        instance->m_has_file_names = true;
    }
    if (function != nullptr) {
        function_mangled = function;
        function_demangled = Backtrace::demangle(function);

        if (function_mangled == "main") { // avoid recording the next frames
            instance->m_skip_next_frames = true;
        }
    }

    // int pc, const std::string& filename, int lineno, const std::string& fn_demangled, const std::string& fn_mangled, const std::string& fn_demangled);
    instance->m_backtrace.push_back(Backtrace::Frame{pc, str_filename, lineno, function_mangled, function_demangled});

    return 0;
}
} // extern "C"

/**
 * Backtrace
 */

Backtrace::Backtrace(size_t frames_to_skip) {
    unique_lock<mutex> lock(g_mutex);

    // State machine, which operation failed/threw the exception
    enum {
        Initialisation,
        Backtrace_full
    } operation;

    try {
        // Initialise the state of the library
        operation = Initialisation;
        if (g_backtrace_state == nullptr) {
            g_backtrace_state = backtrace_create_state(
                    /* program name (NULL => autodeduce) */ nullptr,
                    /* threaded ? */ 0,
                    /* callback in case of failure */ on_error,
                    /* instance */ nullptr);
        }

        // Create the backtrace
        operation = Backtrace_full;
        int rc = backtrace_full(g_backtrace_state,
                                frames_to_skip + 1 /* ignore ctor */, _on_backtrace_full, on_error, this);
        if (rc != 0) { throw runtime_error("backtrace failed"); }

    } catch (std::runtime_error& e) {
        stringstream ss;
        ss << "[ERROR] ";
        switch (operation) {
        case Initialisation:ss << "on initialisation";
            break;
        case Backtrace_full:ss << "on full backtrace";
            break;
        default:ss << "on unknown";
        }
        ss << ": " << e.what();

        m_error = ss.str();
    }
}

const std::vector<Backtrace::Frame>& Backtrace::get_backtrace() const {
    return m_backtrace;
}

const std::string& Backtrace::get_error() const {
    return m_error;
}

const bool Backtrace::has_error() const noexcept {
    return !m_error.empty();
}

const bool Backtrace::has_debug_info() const noexcept {
    return m_has_file_names;
}

/**
* Demangle a given function name
*/
string Backtrace::demangle(const char* fn) {
    int rc = 0;
    string result;
    char* cresult = abi::__cxa_demangle(fn, nullptr, nullptr, &rc);
    if (rc == 0) { // ok
        result = cresult;
    }
    free(cresult);
    cresult = nullptr;

    return result;
}

string Backtrace::demangle(const std::string& fn) {
    return demangle(fn.c_str());
}

/**
 * Backtrace::Entry
 */
Backtrace::Frame::Frame(uint64_t pc, const string& filename, int lineno, const string& fn_mangled,
                        const std::string& fn_demangled) :
        m_program_counter(pc), m_file_name(filename), m_line_number(lineno), m_function_name_mangled(fn_mangled),
        m_function_name_demangled(fn_demangled) { }

uint64_t Backtrace::Frame::get_program_counter() const noexcept {
    return m_program_counter;
}

const string& Backtrace::Frame::get_file_name() const noexcept {
    return m_file_name;
}

int Backtrace::Frame::get_line_number() const noexcept {
    return m_line_number;
}

const string& Backtrace::Frame::get_function_name() const noexcept {
    if (!get_function_name_demangled().empty())
        return get_function_name_demangled();
    else
        return get_function_name_mangled();
}

const string& Backtrace::Frame::get_function_name_mangled() const noexcept {
    return m_function_name_mangled;
}

const string& Backtrace::Frame::get_function_name_demangled() const noexcept {
    return m_function_name_demangled;
}

std::ostream& operator<<(std::ostream& out, const ::common::Backtrace& backtrace) {
    if (backtrace.has_error()) {
        out << backtrace.get_error() << '\n';
    } else {
        auto& bt = backtrace.get_backtrace();
        for (size_t i = 0; i < bt.size(); i++) {
            out << "[#" << i << "] " << bt[i] << '\n';
        }
    }

    return out;
}

std::ostream& operator<<(std::ostream& out, const ::common::Backtrace::Frame& entry){
    out << "pc: 0x" << hex << entry.get_program_counter() << dec << " at " <<  entry.get_file_name() << ":" << entry.get_line_number() << ", function: " <<
        entry.get_function_name();
    return out;
}

} // namespace common