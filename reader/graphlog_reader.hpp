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

#pragma once

#include "reader.hpp"

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <memory>
#include <random>
#include <unordered_map>

namespace gfe::graph { class WeightedEdge; } // forward decl.

namespace gfe::reader::graphlog {

/**
 * A graphlog file consists of multiple sections of content:
 * 1- First is the header, which contains all the metadata and the properties stored by the generator. This is in plain format and it is similar to
 *    a java .properties file with each properties stored as `key = value'. Use the global function graphlog::parse_properties(path_file) to read
 *    the properties from the file
 * 2- There are two compressed sections for the vertices, one is for the final and the other is for the temporary vertices of the graph. Final vertices
 *    are present in the original graph used to create the log, while temporary vertices are artificial entities created by the generator to make some
 *    noise (and gaps) in the update process. At the end of an update, these `temporary' vertices need to be removed.
 *    There are two utility classes to load the vertices from an input file:
 *    a. VertexLoader to load the vertices in `bulk'
 *    b. VertexReader to load one vertex at the time
 *    In both cases, the fstream handle need to be positioned at the start of the compressed section for the final or the temporary vertices. That is
 *    specified in the properties `internal.vertices.final.begin' and `internal.vertices.temporary.begin'. It is not possible to read with the same
 *    instance read both sections (final and temporary).
 * 3- Finally there are `internal.edges.num_blocks' compressed sections consisting of updates. Each section consists of three arrays, one for the sources,
 *    one for the destinations, and one for the weights. A weight has the value -1 if the update represents a deletion, 0 if the weight is not specified
 *    and >0 if the weight is the same of the original graph. Again there are two utility classes to load the edges:
 *    a. EdgeLoader to load the edges in `bulk'. The size of the array must be internal.edges.block_size in bytes and internal.edges.block_size / (3*sizeof(uint64_t))
 *       in terms of num_edges.
 *    b. EdgeReader to read one edge at the time.
 */

// A `graphlog' contains a set of properties, stored in plain format as "name = value". We represent these properties in an hash table.
using Properties = std::unordered_map<std::string, std::string>;

// Read the properties from the given file
Properties parse_properties(const std::string& path_graphlog);

// Read the properties from the given handle
Properties parse_properties(std::fstream& handle);

// Set the current position of the handle at the start of the given section
enum class Section { VTX_FINAL, VTX_TEMP, EDGES };
void set_marker(const Properties& properties, std::fstream& handle, Section section);

// Load the vertices from the given graph. The handle should already be position at the start of the respective compressed section.
class VertexLoader {
    VertexLoader(const VertexLoader&) = delete;
    VertexLoader& operator=(const VertexLoader&) = delete;

    std::fstream& m_handle; // handle to read the vertices from the file
    static constexpr uint64_t m_input_stream_sz = (1ull << 20); // 1 Mb
    uint8_t* m_input_stream { nullptr }; // compressed content read from the handle
    std::streampos m_input_stream_pos { 0 }; // offset of the next chunk to read the file
    void* m_zstream {nullptr}; // current stream for zlib, to read the vertices from the input file
    static constexpr uint64_t m_output_stream_sz = (1ull << 17); // * sizeof(uint64_t)
    uint64_t m_output_leftover { 0 };
    uint64_t m_output_leftover_sz { 0 };


public:
    // Initialise the reader, the handle should already be positioned at the start of the compressed section
    VertexLoader(std::fstream& handle);

    // Destructor
    ~VertexLoader();

    // Load the up to `array_sz' vertices in the given buffer. Return the number of vertices loaded, or 0 if the iterator has been depleted.
    uint64_t load(uint64_t* array, uint64_t array_sz);
};

// Read one vertex at the time from the input file
class VertexReader {
    VertexReader(const VertexLoader&) = delete;
    VertexReader& operator=(const VertexReader&) = delete;

    VertexLoader m_loader;
    static constexpr uint64_t m_buffer_capacity = (1ull << 20); // 1 M vertices, 8 MB
    uint64_t m_buffer_size = 0; // total number of vertices currently in the buffer
    uint64_t* m_buffer { nullptr }; // vertices decompressed by the loader
    uint64_t m_position { 0 }; // next vertex to read from the buffer

public:
    // Initialise the reader, the handle should already be position at the start of the compressed section
    VertexReader(std::fstream& handle);

    // Destructor
    ~VertexReader();

    // Read one vertex at the time
    bool read_vertex(uint64_t& out_vertex);
};


// Loader of whole blocks of edges
class EdgeLoader {
    EdgeLoader(const EdgeLoader&) = delete;
    EdgeLoader& operator=(const EdgeLoader&) = delete;

    std::fstream& m_handle; // handle to read the edges from the file
    static constexpr uint64_t m_input_stream_sz = (1ull << 20); // 256 KB
    uint8_t* m_input_stream { nullptr }; // compressed content read from the handle
    std::streampos m_input_stream_pos { 0 }; // offset of the next chunk to read the file

public:
    // Initialise the reader, the handle should already be positioned at the start of the compressed section
    EdgeLoader(std::fstream& handle);

    // Destructor
    ~EdgeLoader();

    // Load a whole block of edges in the given buffer. Return the number of edges loaded, or 0 if the array is not big enough to load the whole block
    uint64_t load(uint64_t* array, uint64_t num_edges);
};

// Iterate over an edge at the time from a block of edges
class EdgeBlockReader {
    std::shared_ptr<uint64_t[]> m_ptr_block; // the block of edges
    uint64_t m_position; // the current position in the block
    uint64_t m_num_edges; // the total number of edges in the current block

public:
    // Create a dummy reader
    EdgeBlockReader();

    // Create a new instance to iterate over of a block of edges, as retrieved from the EdgeLoader
    EdgeBlockReader(std::shared_ptr<uint64_t[]> block, uint64_t num_edges);

    // Copy ctor and assignment are allowed in this case
    EdgeBlockReader(const EdgeBlockReader&) = default;
    EdgeBlockReader& operator=(const EdgeBlockReader&) = default;

    // Destructor
    ~EdgeBlockReader();

    // Check whether there are more edges to read
    bool has_next() const;

    // Read one edge at the time
    bool read_edge(graph::WeightedEdge& edge);
};

// Read one block of edges at the time
class EdgeBlockLoader {
    EdgeBlockLoader(const EdgeBlockLoader&) = delete;
    EdgeBlockLoader& operator=(const EdgeBlockLoader&) = delete;

    EdgeLoader m_loader;
    const uint64_t m_max_num_edges; // max number of edges that can be stored in the buffer
    std::shared_ptr<uint64_t[]> m_ptr_buffer = 0; // pointer to where the edges are stored

public:
    // Create a new reader. The handle should be already positioned at the start of the compressed stream. The value for max_num_edges
    // should be the same of value of the property `internal.edges.block_size'
    EdgeBlockLoader(std::fstream& handle, uint64_t block_size_bytes);

    // Destructor
    ~EdgeBlockLoader();

    // Load one block of edges from the file
    EdgeBlockReader load();
};

// Read one edge at the time
class EdgeReader : public ::gfe::reader::Reader {
    EdgeReader(const EdgeReader&) = delete;
    EdgeReader& operator=(const EdgeReader&) = delete;

    std::fstream m_handle;
    EdgeBlockLoader m_loader;
    EdgeBlockReader m_reader;

    // Internal ctor
    EdgeReader(const std::string& path, Properties properties);

public:
    // Create a new reader for the given path
    EdgeReader(const std::string& path);

    // Destructor
    virtual ~EdgeReader();

    // Retrieve the next edge of the file. Returns true if an edge has been read, false if we reached the end of the file.
    virtual bool read(graph::WeightedEdge& edge) override;
    bool read_edge(graph::WeightedEdge& edge); // alias

    // All graphlog graphs are currently undirected only
    virtual bool is_directed() const override;
};

} // namespace
