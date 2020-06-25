/*
 * dummy.cpp
 *
 *  Created on: 22 Aug 2019
 *      Author: Dean De Leo
 */

#include "dummy.hpp"

using namespace std;

namespace gfe::library {

Dummy::Dummy(bool is_directed) : m_is_directed(is_directed) { }
Dummy::~Dummy() { }
uint64_t Dummy::num_edges() const { return 0; }
uint64_t Dummy::num_vertices() const { return 0; }
bool Dummy::is_directed() const { return m_is_directed; }
bool Dummy::has_vertex(uint64_t vertex_id) const { return true; }
double Dummy::get_weight(uint64_t source, uint64_t destination) const { return 0.0; }
void Dummy::dump_ostream(std::ostream& out) const { out << "DuMMy"; }
bool Dummy::add_vertex(uint64_t vertex_id) { return true; }
bool Dummy::remove_vertex(uint64_t vertex_id){ return true; }
bool Dummy::add_edge(graph::WeightedEdge e) { return true; }
bool Dummy::add_edge_v2(graph::WeightedEdge e) { return true; }
bool Dummy::remove_edge(graph::Edge e){ return true; }
void Dummy::set_timeout(uint64_t seconds) { }

} // namespace


