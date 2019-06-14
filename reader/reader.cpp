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

#include "dimacs9_reader.hpp"
#include "format.hpp"
#include "metis_reader.hpp"
#include "plain_reader.hpp"

#undef CURRENT_ERROR_TYPE
#define CURRENT_ERROR_TYPE reader::ReaderError

using namespace std;

namespace reader {

Reader::Reader(){ }
Reader::~Reader(){ }

std::unique_ptr<Reader> Reader::open(const std::string& path){
    auto format = get_graph_format(path);
    switch(format){
    case Format::DIMACS9:
        return make_unique<Dimacs9Reader>(path);
    case Format::METIS:
        return make_unique<MetisReader>(path);
    case Format::PLAIN:
        return make_unique<PlainReader>(path, false);
    case Format::PLAIN_WEIGHTED:
        return make_unique<PlainReader>(path, true);
    default:
        ERROR("Unrecognised graph format for the file: `" << path << "'");
    }
}

} // namespace reader

