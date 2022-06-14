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

#include "error.hpp"

#include <iostream>

namespace common {

/**
 * Remove the starting colons (::) from the exception class
 */
static std::string trim_exception_class(const std::string &class_name_arg) {
    std::string class_name = class_name_arg;
    if (class_name.find("::") == 0) { // remove the prefix :: for the global namespace
        class_name = class_name.substr(2);
    }
    if (class_name.find("common::") == 0) { // remove the prefix (namespace) common::
        class_name = class_name.substr(8);
    }
    return class_name;
}


Error::Error(const std::string &exceptionClass_, const std::string &message_, const std::string &file_, int line_, const std::string &function_)
        : runtime_error(message_), m_class(trim_exception_class(exceptionClass_)), m_file(file_), m_line_number(line_), m_function_name(function_),
          m_backtrace(/* skip Error ctor */ 1) {}

std::string Error::get_file() const { return m_file; }

int Error::get_line_number() const { return m_line_number; }

std::string Error::get_function_name() const { return m_function_name; }

std::string Error::get_exception_class() const { return m_class; }

const Backtrace &Error::get_backtrace() const { return m_backtrace; }

// Definition of the utility stream
thread_local std::stringstream Error::utilitystream;

std::ostream &operator<<(std::ostream &out, const Error &e) {
    out << "Kind: " << e.get_exception_class() << ", file: " << e.get_file() << ", function: " << e.get_function_name()
        << ", line: " << e.get_line_number() << "\n";
    auto &bt = e.get_backtrace();
    if (!bt.has_debug_info()) {
        out << "> Backtrace not available as debug symbols (-g) are not present.\n";
    } else {
        out << "Backtrace:\n" << bt;
    }
    out << "ERROR: " << e.what();
    return out;
}


} // namespace common
