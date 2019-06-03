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

#include "plain_reader.hpp"

#include <fstream>
#include <string>
#include "graph/edge.hpp"
#include "configuration.hpp"
#include "utility.hpp"

using namespace std;

#undef CURRENT_ERROR_TYPE
#define CURRENT_ERROR_TYPE reader::ReaderError

namespace reader {

PlainWeightedReader::PlainWeightedReader(const string& path) : m_handle(init_fstream(path)){
    LOG("[PlainWeightedReader] Reading `" << path << "' ...");
}

PlainWeightedReader::~PlainWeightedReader(){
    m_handle.close();
}

bool PlainWeightedReader::read(graph::WeightedEdge& edge) {
    if(m_handle.eof()) return false;
    m_handle >> edge.m_source >> edge.m_destination >> edge.m_weight;
    if(m_handle.fail()) ERROR("Error while processing the input file");
    return true;
}

} // namespace reader
