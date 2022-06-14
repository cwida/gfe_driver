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

#ifndef COMMON_BACKTRACE_HPP
#define COMMON_BACKTRACE_HPP

#include <cinttypes>
#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

namespace common {

// forward declaration, implementation detail
extern "C" {
int _on_backtrace_full(void* instance, uintptr_t pc, const char* filename, int lineno, const char* function);
};

/**
 *  Record a back trace from the given program counter. It mostly works when the build includes the debug symbols,
 *  otherwise the function names are not recognised.
 */
class Backtrace {
    friend int
    ::common::_on_backtrace_full(void* instance, uintptr_t pc, const char* filename, int lineno, const char* function);

public:
    /**
     * A single entry in the back trace.
     */
    class Frame {
        friend class Backtrace;

        uint64_t m_program_counter;
        std::string m_file_name;
        int m_line_number;
        std::string m_function_name_mangled;
        std::string m_function_name_demangled;
    public:

        /**
         * Initialise the entry
         * @param pc current program counter
         * @param filename full path to the C++ source file where the execution relates
         * @param lineno position in `filename' for the current execution
         * @param fn_mangled the function name related to the program counter. For C++ the function name is mangled
         * @param fn_demangled the function name related to the program counter. For C++ this is the full function signature, demangled.
         */
        Frame(uint64_t pc, const std::string& filename, int lineno, const std::string& fn_mangled,
              const std::string& fn_demangled);

        /**
         * Retrieve the registered program counter
         */
        uint64_t get_program_counter() const noexcept;

        /**
         * Retrieve the registered file name, if present
         */
        const std::string& get_file_name() const noexcept;

        /**
         * Retrieve the line number associated to the file name, where the frame points to.
         */
        int get_line_number() const noexcept;

        /**
         * Retrieve the name of the function where the program counter is. It returns the demangled signature if available,
         * otherwise the mangled one.
         */
        const std::string& get_function_name() const noexcept;

        /**
         * Retrieve the function name, in its original form (mangled).
         */
        const std::string& get_function_name_mangled() const noexcept;

        /**
         * Retrieve the full signature demangled, if present.
         */
        const std::string& get_function_name_demangled() const noexcept;

    };

private:
    std::vector<Frame> m_backtrace; // the list of frames registered
    std::string m_error; // in case an error has occurred,
    bool m_has_file_names = false; // whether the backtrace recorded any file name. Typically true for Debug builds (-g).
    bool m_skip_next_frames = false; // whether to avoid recording the next frames during a backtrace


public:
    /**
     * Creates a new back trace from the current frame
     * @param frames_to_skip the number of top frame to ignore in the recording
     */
    Backtrace(size_t frames_to_skip = 0);

    /**
     * Retrieve the registered backtrace
     */
    const std::vector<Frame>& get_backtrace() const;

    /**
     * Checks whether the current backtrace is in an erroneous state and a recording has not been properly registered.
     */
    const bool has_error() const noexcept;

    /**
     * In case of error, it returns a message with the cause of the error.
     */
    const std::string& get_error() const;


    /**
     * Checks whether the debug symbols are present. Without debug symbols a backtrace only contains a list
     * of program counters for the frames associated to the execution.
     */
    const bool has_debug_info() const noexcept;


    /**
     * Demangle the given function name
     */
    static std::string demangle(const char* fn);

    static std::string demangle(const std::string& fn);

};

/**
 * Output
 */
std::ostream& operator<<(std::ostream& out, const common::Backtrace::Frame& entry);
std::ostream& operator<<(std::ostream& out, const common::Backtrace& backtrace);

} // namespace common
#endif // COMMON_BACKTRACE_HPP
