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

#ifndef COMMON_ERROR_HPP
#define COMMON_ERROR_HPP

#include <ostream>
#include <sstream>
#include <stdexcept>

#include "backtrace.hpp"

namespace common {

/**
 * Base class for any errors
 */
class Error : public std::runtime_error {
private:
    const std::string m_class; // the class of the exception (eg. "Graph::Exception")
    const std::string m_file; // source file where the exception has been raised
    int m_line_number; // related line where the exception has been raised
    const std::string m_function_name; // function causing the exception
    Backtrace m_backtrace;

public:
    /**
     * Constructor
     * @param exceptionClass the name of the exception class
     * @param message the error message associated to this exception
     * @param file the source file where the exception has been generated
     * @param line the line where the exception has been generated
     * @param function the function where this exception has been raised
     */
    Error(const std::string& exceptionClass, const std::string& message, const std::string& file, int line,
          const std::string& function);

    /**
     * Retrieve the source file where the exception has been raised
     * @return the source file where the exception has been raised
     */
    std::string get_file() const;

    /**
     * The line number where the exception has been raised
     * @return the line number where the exception has been raised
     */
    int get_line_number() const;

    /**
     * Retrieve the function that fired the exception
     * @return the function name where the exception has been raised
     */
    std::string get_function_name() const;

    /**
     * Retrieves the name of the exception class
     * @return the name of the exception class
     */
    std::string get_exception_class() const;

    /**
     * Retrieve the backtrace for the given exception.
     */
    const Backtrace& get_backtrace() const;

    /**
     * Utility class to create the exception message
     */
    static thread_local std::stringstream utilitystream;
};

// In case of changing the inheritance of the exception and the macro is caught from somewhere else
#if !defined(RAISE_EXCEPTION_CREATE_ARGUMENTS)
/**
 * It prepares the arguments `file', `line', `function', and `what' to be passed to an exception ctor
 * @param msg the message stream to concatenate
 */
#define RAISE_EXCEPTION_CREATE_ARGUMENTS(msg) const char* file = __FILE__; int line = __LINE__; const char* function = __FUNCTION__; \
        auto& stream = ::common::Error::utilitystream; \
        stream.str(""); stream.clear(); \
        stream << msg; \
        std::string what = stream.str(); \
        stream.str(""); stream.clear() /* reset once again */
#endif

/**
 * Raises an exception with the given message
 * @param exception the exception to throw
 * @param msg: an implicit ostream, with arguments concatenated with the symbol <<
 */

#define RAISE_EXCEPTION(exc, msg) { RAISE_EXCEPTION_CREATE_ARGUMENTS(msg); throw exc( #exc, what, file, line, function); }


/**
 * These exception classes are so similar, so define a general macro to create the exception
 */
#define DEFINE_EXCEPTION_SUBCLASS(subClass, superClass) class subClass: public superClass { \
        public: subClass(const std::string& exceptionClass, const std::string& message, const std::string& file, \
            int line, const std::string& function) : \
                superClass(exceptionClass, message, file, line, function) { } \
} /* End of DEFINE_EXCEPTION_SUBCLASS */
#define DEFINE_EXCEPTION(exceptionName) DEFINE_EXCEPTION_SUBCLASS(exceptionName, ::common::Error)

/**
 * The exception type to throw when using the macro `ERROR'
 */
#define CURRENT_ERROR_TYPE ::common::Error

/**
 * Helper for the macro `ERROR' to fully expand the definition CURRENT_ERROR_TYPE before
 * invoking RAISE_EXCEPTION
 */
#define _RAISE_EXCEPTION(exc, msg) RAISE_EXCEPTION(exc, msg)

/**
 * Convenience macro, it raises an exception of type `Error'.
 */
#define ERROR(message) _RAISE_EXCEPTION(CURRENT_ERROR_TYPE, message)

/**
 * Convenience macro, runtime assertion
 */
#define ASSERT(condition) if(!(condition)){ ERROR("Assertion failed: " << #condition); }


/**
 * Standard exception when an invalid argument has been given to a function, similar to std::invalid_argument
 */
DEFINE_EXCEPTION(InvalidArgument);
#define INVALID_ARGUMENT(message) RAISE_EXCEPTION(::common::InvalidArgument, message)

/**
 * Overload the operator to print the descriptive content of an ELF Exception
 */
std::ostream& operator<<(std::ostream& out, const ::common::Error& e);

} // namespace common


#endif /* COMMON_ERROR_HPP */
