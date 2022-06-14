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


#ifndef COMMON_FILESYSTEM_HPP
#define COMMON_FILESYSTEM_HPP

#include <string>

namespace common::filesystem {

/**
 * Change the current w.d. to the specified path, and restore the previous
 * working directory on exit.
 */
class TemporaryWorkingDirectory{
private:
    const std::string m_old_wd; // the working directory to restore on exit
public:

    /**
     * Restore the current working directory on exit, but do not explicitly change it
     */
    TemporaryWorkingDirectory();

    /**
     * Change the working directory to the given path
     * @param path path to the new working directory
     */
    TemporaryWorkingDirectory(const std::string& path);

    /**
     * Restore the previous working directory
     */
    ~TemporaryWorkingDirectory() noexcept(false);
};


/**
 * Get the absolute path for the given (relative) path
 */
std::string absolute_path(const std::string& path);


/**
 * Get the directory path for the given path
 */
std::string directory(const std::string& path);

/**
 * Check whether the given file exists
 */
bool exists(const std::string& path);
bool file_exists(const std::string& path);

/**
 * Check whether the given path exists and it is a directory
 */
bool is_directory(const std::string& path);

/**
 * Retrieve the size, in bytes, of the file in the given path
 */
uint64_t file_size(const std::string& path);

/**
 * Retrieve the base name of the given path
 */
std::string filename(const std::string& path);

/**
 * Retrieve the extension of the given path
 */
std::string extension(const std::string& path);

/**
 * Get the path to the current working directory
 */
std::string wd();
std::string working_directory();

/**
 * Retrieve the absolute path to the program executable
 */
std::string path_executable();

/**
 * Retrieve the absolute path to the directory of the executable
 */
std::string directory_executable();

/**
 * Create the given directory. All the parents are also created in the path,
 * emulating the same functionality of `mkdir -pv'
 */
void mkdir(const std::string& path);

} // common::filesystem

#endif //COMMON_FILESYSTEM_HPP
