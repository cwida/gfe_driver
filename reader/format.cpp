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

#include "format.hpp"

#include <cstring>
#include <strings.h>
#include "common/filesystem.hpp"

namespace reader {

Format get_graph_format(const char* path) {
    const char* extension = strrchr(path, '.');
    if(extension != nullptr){
        extension++; // skip the '.'
        if ( strcasecmp(extension, "properties") == 0 ){ // properties file of the Graphalytics suite
            return reader::Format::LDBC_GRAPHALYTICS;
        } else if( strcasecmp(extension, "log") == 0 || strcasecmp(extension, "glog") == 0 || strcasecmp(extension, "graphlog") == 0){
            return reader::Format::GRAPHLOG;
        } else if( strcasecmp(extension, "el") == 0){
            return reader::Format::PLAIN;
        } else if( strcasecmp(extension, "wel") == 0 ){
            return reader::Format::PLAIN_WEIGHTED;
        } else if( strcasecmp(extension, "metis") == 0 || strcasecmp(extension, "graph") == 0) {
            return reader::Format::METIS;
        } else if( strcasecmp(extension, "gr") == 0 || strcasecmp(extension, "dimacs") == 0 || strcasecmp(extension, "dimacs9") == 0 ){
            return reader::Format::DIMACS9;
        }
    }

    return reader::Format::UNKNOWN;
}

Format get_graph_format(const std::string& path) {
    return get_graph_format(path.c_str());
}


} // namespace reader
