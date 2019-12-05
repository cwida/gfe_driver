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

#include "reader.hpp"
#include "utility.hpp"

#include <fstream>
#include "common/filesystem.hpp"

using namespace std;

#undef CURRENT_ERROR_TYPE
#define CURRENT_ERROR_TYPE ::gfe::reader::ReaderError

namespace gfe::reader {

fstream init_fstream(const string& path){
    fstream handle(path.c_str(), ios_base::in);
    if(!handle.good()){ // some error occurred
        if(!common::filesystem::file_exists(path)){
            ERROR("The file `" << path << "' does not exist");
        } else {
            ERROR("Cannot read the file: `" << path << "'");
        }
    }
    return handle;
}

} // namespace reader
