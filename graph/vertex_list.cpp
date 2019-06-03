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

#include "vertex_list.hpp"

#include "common/error.hpp"
#include "common/permutation.hpp"
#include "common/timer.hpp"
#include "cbytearray.hpp"
#include "configuration.hpp"

using namespace common;
using namespace std;

namespace graph {

VertexList::VertexList(CByteArray* vertices) : m_vertices(vertices) {
    ASSERT(vertices != nullptr && "The given argument is a nullptr");
}
VertexList::~VertexList(){ delete m_vertices; m_vertices = nullptr; }

uint64_t VertexList::get(uint64_t index) const {
    if(index >= num_vertices()) INVALID_ARGUMENT("The given index is out of bounds: " << index << " >= " << num_vertices());
    return m_vertices->get_value_at(index);
}

uint64_t VertexList::num_vertices() const noexcept {
    return m_vertices->capacity();
}

void VertexList::permute() {
    permute(configuration().seed() + 106);
}

void VertexList::permute(uint64_t seed){
    LOG("Permuting the list of " << num_vertices() << " vertices ...");
    Timer timer; timer.start();

    // create the permutation array
    auto ptr_permutation = make_unique<uint64_t[]>(num_vertices());
    uint64_t* __restrict permutation = ptr_permutation.get();
    for(uint64_t i = 0, N = num_vertices(); i < N; i++) permutation[i] = i;
    common::permute(permutation, num_vertices(), seed);

    auto new_vertices = make_unique<CByteArray>(m_vertices->get_bytes_per_element(), num_vertices());
    for(uint64_t i = 0, N = num_vertices(); i < N; i++){
        new_vertices->set_value_at(i, m_vertices->get_value_at(permutation[i]));
    }

    timer.stop();
    LOG("List of vertices permuted in " << timer);

    delete m_vertices; m_vertices = new_vertices.release();
}

} // namespace graph
