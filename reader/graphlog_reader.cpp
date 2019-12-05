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

#include "graphlog_reader.hpp"

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <regex>
#include <stdexcept>
#include "zlib.h"

#include "common/filesystem.hpp"
#include "graph/edge.hpp"
#include "configuration.hpp"

using namespace std;

#undef CURRENT_ERROR_TYPE
#define CURRENT_ERROR_TYPE ::gfe::reader::ReaderError

namespace gfe::reader {

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG
#define COUT_DEBUG_FORCE(msg) { std::cout << "[Graphlog::" << __FUNCTION__ << "] " << msg << std::endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif

/*****************************************************************************
 *                                                                           *
 *  Properties                                                               *
 *                                                                           *
 *****************************************************************************/
namespace graphlog {

Properties parse_properties(const std::string& path_graphlog){
    if(!common::filesystem::file_exists(path_graphlog)) ERROR("The given file does not exist: " << path_graphlog);
    fstream handle{path_graphlog, ios_base::in | ios_base::binary};
    if(!handle.good()) ERROR("Cannot open the file: " << path_graphlog);
    Properties properties = parse_properties(handle);
    handle.close();

    return properties;
}

// Read the properties from the given handle
Properties parse_properties(std::fstream& handle){
    Properties properties;
    string line;
    getline(handle, line);
    if(line != "# GRAPHLOG") ERROR("Missing magic header '# GRAPHLOG'");

    regex pattern { "^\\s*([A-Za-z0-9_.-]+?)\\s*=\\s*([^#\n]+?)\\s*" };
    while(!handle.eof()){
        string line;
        getline(handle, line);
        if(line == "__BINARY_SECTION_FOLLOWS") break; // done

        smatch matches;
        if ( regex_match(line, matches, pattern) ) { // side effect => populate matches
            string key = matches[1];
            string value = matches[2];

            properties[key] = value;
        }
    }

    return properties;
}

void set_marker(const Properties& properties, std::fstream& handle, Section section){
    Properties::const_iterator property;

    switch(section){
    case Section::VTX_FINAL:
        property = properties.find("internal.vertices.final.begin");
        break;
    case Section::VTX_TEMP:
        property = properties.find("internal.vertices.temporary.begin");
        break;
    case Section::EDGES:
        property = properties.find("internal.edges.begin");
        break;
    }

    if(property == properties.end()) ERROR("Missing required property");

    handle.clear(); // if eof() -> clear the flags
    handle.seekg(stoull(property->second));
}

} // namespace graphlog

/*****************************************************************************
 *                                                                           *
 *  VertexLoader                                                             *
 *                                                                           *
 *****************************************************************************/
namespace graphlog {

VertexLoader::VertexLoader(std::fstream& handle) : m_handle(handle) {
    m_input_stream_pos = m_handle.tellg();

    // internal buffers
    m_input_stream = new uint8_t[m_input_stream_sz];

    // initialise the compressed stream
    m_zstream = malloc(sizeof(z_stream));
    if(m_zstream == nullptr) throw bad_alloc();
    z_stream* zstream = reinterpret_cast<z_stream*>(m_zstream);
    zstream->zalloc = Z_NULL;
    zstream->zfree = Z_NULL;
    zstream->opaque = Z_NULL;
    zstream->avail_in = 0;
    zstream->next_in = (unsigned char*) m_input_stream;
    zstream->avail_out = 0;
    zstream->next_out = nullptr;
    int rc = inflateInit2(zstream, -15);
    if(rc != Z_OK) ERROR("Cannot initialise the library zlib");
}

VertexLoader::~VertexLoader(){
    if(m_zstream != nullptr){
        inflateEnd(reinterpret_cast<z_stream*>(m_zstream));
        free(m_zstream); m_zstream = nullptr;
    }

    delete[] m_input_stream; m_input_stream = nullptr;
}

uint64_t VertexLoader::load(uint64_t* array, uint64_t array_sz){
    if(array == nullptr) INVALID_ARGUMENT("The argument `array' is null");
    if(m_zstream == nullptr || array_sz == 0) return 0;
    z_stream* zstream = reinterpret_cast<z_stream*>(m_zstream);
    uint64_t num_elements_loaded = 0;
    bool depleted = false;

    do {
        // read the content from the input file
        if(zstream->avail_in == 0){
            m_handle.seekp(m_input_stream_pos);
            m_handle.read((char*) m_input_stream, m_input_stream_sz);
            uint64_t data_read = m_input_stream_sz;
            if(m_handle.eof()){
                data_read = static_cast<uint64_t>(m_handle.tellg()) - m_input_stream_sz;
            } else if(m_handle.bad()){ ERROR("Cannot read from the input file"); }

            zstream->next_in = m_input_stream;
            zstream->avail_in = data_read;
            m_input_stream_pos += data_read;

            COUT_DEBUG("[INPUT VERTICES] " << (int) m_input_stream[0] << ":" << (int) m_input_stream[1] << ":" << (int) m_input_stream[2] << ":" << (int) m_input_stream[3]);
        }

        // copy the remaining m_output_leftover at the start of the stream
        array[0] = m_output_leftover;
        uint8_t* output_buffer = (uint8_t*) array;
        uint64_t output_buffer_sz = array_sz * sizeof(uint64_t) - m_output_leftover_sz;

        // decompress the input
        zstream->next_out = output_buffer + m_output_leftover_sz;
        zstream->avail_out = output_buffer_sz;
        int rc = inflate(zstream, Z_NO_FLUSH);
        if(rc != Z_OK && rc != Z_STREAM_END) ERROR("Cannot decompress the input stream: " << zstream->msg << " (rc: " << rc << ")");
        depleted = (rc == Z_STREAM_END);
        uint64_t data_processed = m_output_leftover_sz + (output_buffer_sz - zstream->avail_out);
        num_elements_loaded = data_processed / sizeof(uint64_t);
        m_output_leftover_sz = data_processed % sizeof(uint64_t);
        if(m_output_leftover_sz > 0) { m_output_leftover = array[num_elements_loaded]; }
        assert(depleted == false || m_output_leftover_sz == 0);
        COUT_DEBUG("data_processed: " << data_processed << ", first value: " << array[0]);
    } while (num_elements_loaded == 0 && !depleted);

    if(depleted){
        inflateEnd(zstream);
        free(m_zstream); m_zstream = nullptr;
    }

    return num_elements_loaded;
}

} // namespace

/*****************************************************************************
 *                                                                           *
 *  VertexReader                                                             *
 *                                                                           *
 *****************************************************************************/

namespace graphlog {

VertexReader::VertexReader(fstream& handle) : m_loader(handle), m_buffer(new uint64_t[m_buffer_capacity]) {

}

VertexReader::~VertexReader(){
    delete[] m_buffer; m_buffer = nullptr;
}

bool VertexReader::read_vertex(uint64_t& out_vertex){
    if(m_position >= m_buffer_size){
        m_buffer_size = m_loader.load(m_buffer, m_buffer_capacity);
        if(m_buffer_size == 0) return false;
    }

    out_vertex = m_buffer[m_position];
    m_position++;
    return true;
}

} // namespace

/*****************************************************************************
 *                                                                           *
 *  EdgeLoader                                                               *
 *                                                                           *
 *****************************************************************************/

namespace graphlog {
EdgeLoader::EdgeLoader(std::fstream& handle) : m_handle(handle) {
    m_input_stream = new uint8_t[m_input_stream_sz];
}

EdgeLoader::~EdgeLoader(){
    delete[] m_input_stream; m_input_stream = nullptr;
}

uint64_t EdgeLoader::load(uint64_t* array, uint64_t num_edges){
    if(!m_handle.good()) return 0;

    std::streampos handle_pos_start = m_handle.tellg();
    uint64_t num_edges_loaded = 0;
    bool done = false;

    // init the zlib stream
    z_stream z;
    z.zalloc = Z_NULL;
    z.zfree = Z_NULL;
    z.opaque = Z_NULL;
    z.avail_in = 0;
    z.next_in = (unsigned char*) m_input_stream;
    uint64_t output_buffer_sz = num_edges * sizeof(uint64_t) * 3;
    z.avail_out = output_buffer_sz;
    z.next_out = (unsigned char*) array;
    int rc = inflateInit2(&z, -15);
    if(rc != Z_OK) ERROR("Cannot initialise the library zlib");

    do {
        // the input buffer must have been consumed
        assert(z.avail_in == 0);
        uint64_t input_stream_sz = 0;

        // read the content from the input file
        std::streampos handle_pos_last = m_handle.tellg();

        m_handle.read((char*) m_input_stream, m_input_stream_sz);
        input_stream_sz = m_input_stream_sz;
        if(m_handle.eof()){
            m_handle.clear();
            m_handle.seekg(0, ios_base::end);
            input_stream_sz = static_cast<uint64_t>(m_handle.tellg()) - handle_pos_last;
        } else if(m_handle.bad()){ ERROR("Cannot read from the input file"); }

        if(input_stream_sz == 0) break; // EOF, there is nothing to read

        z.next_in = m_input_stream;
        z.avail_in = input_stream_sz;

        COUT_DEBUG("[INPUT EDGES] " << (int) m_input_stream[0] << ":" << (int) m_input_stream[1] << ":" << (int) m_input_stream[2] << ":" << (int) m_input_stream[3] << " z.avail_in: " << z.avail_in);

        // decompress the input
        rc = inflate(&z, Z_NO_FLUSH);
        if(rc != Z_OK && rc != Z_STREAM_END) ERROR("Cannot decompress the input stream: rc: " << rc << ")");

        // there is not enough space to load the whole block in the buffer
        done = z.avail_out == 0 || rc == Z_STREAM_END;

        if(z.avail_out == 0 && rc != Z_STREAM_END){
            m_handle.seekg(handle_pos_start);
        } else if (rc == Z_STREAM_END){
            assert((output_buffer_sz - z.avail_out) % (3 * sizeof(uint64_t)) == 0);
            num_edges_loaded = (output_buffer_sz - z.avail_out) / (3 * sizeof(uint64_t));
            m_handle.seekg(static_cast<int64_t>(handle_pos_last) + (input_stream_sz - z.avail_in), ios_base::beg);
        }
    } while (!done);

    inflateEnd(&z);
    return num_edges_loaded;
}

}

/*****************************************************************************
 *                                                                           *
 *  EdgeBlockReader                                                          *
 *                                                                           *
 *****************************************************************************/

namespace graphlog {

EdgeBlockReader::EdgeBlockReader() : m_ptr_block(), m_position(0), m_num_edges(0) {

}

EdgeBlockReader::EdgeBlockReader(shared_ptr<uint64_t[]> block, uint64_t num_edges) : m_ptr_block(block), m_position(0), m_num_edges(num_edges){

}

EdgeBlockReader::~EdgeBlockReader(){
    /* nop */
}

bool EdgeBlockReader::read_edge(graph::WeightedEdge& edge){
    if(m_position >= m_num_edges){
        return false;
    } else {
        uint64_t* __restrict sources = m_ptr_block.get();
        uint64_t* __restrict destinations = sources + m_num_edges;
        double* __restrict weights = reinterpret_cast<double*>(destinations + m_num_edges);

        edge.m_source = sources[m_position];
        edge.m_destination = destinations[m_position];
        edge.m_weight = weights[m_position];

        m_position++;
        return true;
    }
}

bool EdgeBlockReader::has_next() const {
    return m_position < m_num_edges;
}

} // namespace

/*****************************************************************************
 *                                                                           *
 *  EdgeBlockLoader                                                          *
 *                                                                           *
 *****************************************************************************/
namespace graphlog {

EdgeBlockLoader::EdgeBlockLoader(std::fstream& handle, uint64_t block_size_bytes) : m_loader(handle), m_max_num_edges(block_size_bytes / (3* sizeof(uint64_t))) {
    if(block_size_bytes % (3*sizeof(uint64_t)) != 0) INVALID_ARGUMENT("Invalid block size: " << block_size_bytes);
    m_ptr_buffer.reset(new uint64_t[3 * m_max_num_edges]);
}

EdgeBlockLoader::~EdgeBlockLoader(){
    /* nop */
}

EdgeBlockReader EdgeBlockLoader::load(){
    uint64_t num_edges = m_loader.load(m_ptr_buffer.get(), m_max_num_edges);
    return EdgeBlockReader{m_ptr_buffer, num_edges};
}

} // namespace

/*****************************************************************************
 *                                                                           *
 *  EdgeReader                                                               *
 *                                                                           *
 *****************************************************************************/

namespace graphlog {

EdgeReader::EdgeReader(const std::string& path) : EdgeReader(path, parse_properties(path)) {

}

EdgeReader::EdgeReader(const std::string& path, Properties properties) : m_loader(m_handle, stoull(properties["internal.edges.block_size"])) {
    m_handle.open(path, ios_base::in | ios_base::binary);
    set_marker(properties, m_handle, Section::EDGES);
    m_reader = m_loader.load();
}

EdgeReader::~EdgeReader(){
    m_handle.close();
}

bool EdgeReader::read(graph::WeightedEdge& edge){
    return read_edge(edge);
}

bool EdgeReader::read_edge(graph::WeightedEdge& edge){
    if(!m_reader.has_next()){ // retrieve the next block of edges
        m_reader = m_loader.load();
    }

    return m_reader.read_edge(edge);
}

bool EdgeReader::is_directed() const {
    return false;
}

} // namespace

} // namespace
