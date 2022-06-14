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

#include "graphalytics_reader.hpp"

#include <cassert>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <regex>
#include <zlib.h>
#include "common/filesystem.hpp"
#include "graph/edge.hpp"
#include "configuration.hpp"
#include "utility.hpp"

using namespace std;
using namespace gfe::reader::details;

#undef CURRENT_ERROR_TYPE
#define CURRENT_ERROR_TYPE ::gfe::reader::ReaderError

namespace gfe::reader {

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
//#define DEBUG
#define COUT_DEBUG_FORCE(msg) { std::cout << "[GraphalyticsReader::" << __FUNCTION__ << "] " << msg << std::endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif

/*****************************************************************************
 *                                                                           *
 *  Utility functions                                                        *
 *                                                                           *
 *****************************************************************************/

// Check whether the given line is a comment or empty, that is, it starts with a sharp symbol # or contains no symbols
static bool ignore_line(const char* line){
    if(line == nullptr) return true;
    while(isspace(line[0])) line++;
    return line[0] == '#' || line[0] == '\0';
}
static bool ignore_line(const std::string& line){
    return ignore_line(line.c_str());
}

// Check whether the current marker points to a number
static bool is_number(const char* marker){
    return marker != nullptr && (marker[0] >= '0' && marker[0] <= '9');
}

/*****************************************************************************
 *                                                                           *
 *  Reader implementations                                                   *
 *                                                                           *
 *****************************************************************************/
namespace details {

// base class, interface
class GraphalyticsReaderBaseImpl{
public:
    virtual ~GraphalyticsReaderBaseImpl(){ }

    // Read one edge at the time
    virtual bool read_edge(uint64_t& source, uint64_t& destination, double& weight) = 0;

    // Read one vertex at the time
    virtual bool read_vertex(uint64_t& vertex_id) = 0;
};

// text format, i.e. the original format mandated by the Graphalytics specification
class GraphalyticsPlainReader : public GraphalyticsReaderBaseImpl {
private:
    fstream m_handle_vertex_file; // I/O handle to parse the vertex file
    fstream m_handle_edge_file; // I/O handle to parse the edge file
    const bool m_is_weighted; // whether the edge list contains weights

public:
    GraphalyticsPlainReader(const string& path_vertex_file, const string& path_edge_file,bool is_weighted) : m_is_weighted(is_weighted) {
        COUT_DEBUG("Opening the vertex file `" << path_vertex_file << "'");
        m_handle_vertex_file.open(path_vertex_file, ios::in);

        COUT_DEBUG("Opening the edge file `" << path_edge_file << "'");
        m_handle_edge_file.open(path_edge_file, ios::in);
    }

    virtual ~GraphalyticsPlainReader(){
        m_handle_edge_file.close();
        m_handle_vertex_file.close();
    }

    // read one edge at the time from the edge file
    bool read_edge(uint64_t& source, uint64_t& destination, double& weight) override {
        if(!m_handle_edge_file.good()) return false;

        // read the next line that is not a comment
        string current_line;
        bool skip { true } ;
        do {
            getline(m_handle_edge_file, current_line);

            skip = ignore_line(current_line);
#if defined(DEBUG)
            if(skip) { COUT_DEBUG("line: `" << current_line << "' is a comment or an empty line, skipped"); }
#endif
        } while (skip && m_handle_edge_file.good());
        if(skip) return false; // ended with a comment
        COUT_DEBUG("Parse line: `" << current_line << "'");

        // read the source
        char* next { nullptr };
        const char* current = current_line.c_str();
        while(isspace(current[0])) current++;
        if(!is_number(current)) ERROR("line: `" << current_line << "', cannot read the source vertex");
        source = strtoull(current, &next, /* base */ 10);

        while(isspace(next[0])) next++;
        current = next;
        if(!is_number(current)) ERROR("line: `" << current_line << "', cannot read the destination vertex");
        destination = strtoull(current, &next, 10);

        if(m_is_weighted){
            while(isspace(next[0])) next++;
            current = next;
            if(!is_number(current)) ERROR("line: `" << current_line << "', cannot read the weight");
            weight = strtod(current, nullptr);
        } else {
            weight = 0;
        }
        return true;
    }

    // read one vertex at the time from the vertex file
    bool read_vertex(uint64_t& out_vertex) override {
        out_vertex = 0; // init

        if(!m_handle_vertex_file.good()) return false;

        // read the next line that is not a comment
        string current_line;
        bool skip { true } ;
        do {
            getline(m_handle_vertex_file, current_line);

            skip = ignore_line(current_line);
#if defined(DEBUG)
            if(skip) { COUT_DEBUG("line: `" << current_line << "' is a comment or an empty line, skipped"); }
#endif
        } while (skip && m_handle_vertex_file.good());
        if(skip) return false; // ended with a comment
        COUT_DEBUG("Parse line: `" << current_line << "'");

        const char* current = current_line.c_str();
        while(isspace(current[0])) current++;
        if(!is_number(current)) ERROR("line: `" << current_line << "', cannot read the vertex id");
        out_vertex = strtoull(current, nullptr, /* base */ 10);

        return true;
    }

};

template<typename T>
class GraphalyticsZlibDecompressInput {
    fstream m_handle; // I/O handle to read the data from the input file
    z_stream m_stream; // zlib handle to decompress the content
    constexpr static uint64_t m_buffer_capacity = (1<<15) * sizeof(T); // the capacity of the handles, in bytes
    uint8_t* m_input_buffer { nullptr }; // the content read from the input vertex file
    uint8_t* m_output_buffer { nullptr }; // the content decompressed from the library
    uint64_t m_output_pos = 0; // current position in the output buffer
    uint64_t m_output_sz = 0; // number of vertices in the output buffer
    uint64_t m_output_leftover_sz = 0; // bytes decompressed from the output that do not reach sizeof(T)

    template<bool bogus = true>
    typename std::enable_if<(bogus && sizeof(T) <= 8), void>::type assign(T* __restrict to, T* __restrict from){
        *to = *from;
    }

    template<bool bogus = true>
    typename std::enable_if<(bogus && sizeof(T) > 8), void>::type assign(T* to, T* from){
        memcpy(to, from, sizeof(T));
    }

public:
    GraphalyticsZlibDecompressInput(const string& path_vertex_file){
        m_handle.open(path_vertex_file, ios::in | ios::binary);
        if(!m_handle.good()){ ERROR("Cannot open the input file: " << path_vertex_file); }

        m_input_buffer = new uint8_t[m_buffer_capacity];
        m_output_buffer = new uint8_t[m_buffer_capacity];

        m_stream.zalloc = Z_NULL;
        m_stream.zfree = Z_NULL;
        m_stream.opaque = Z_NULL;
        m_stream.avail_in = 0;
        m_stream.next_in = m_input_buffer;
        m_stream.avail_out = 0;
        m_stream.next_out = nullptr;
        int rc = inflateInit(&m_stream);
        if(rc != Z_OK) ERROR("Cannot initialise the library zlib");
    }

    ~GraphalyticsZlibDecompressInput(){
        inflateEnd(&m_stream); // ignore rc
        m_handle.close();
        delete[] m_input_buffer; m_input_buffer = nullptr;
        delete[] m_output_buffer; m_output_buffer = nullptr;
    }

    bool read(T* item){
        if(m_output_pos >= m_output_sz){
            // read the content from the input file
            if(m_stream.avail_in == 0){
                if(!m_handle.good()) return false; // depleted

                m_handle.read((char*) m_input_buffer, m_buffer_capacity);
                uint64_t data_read = m_handle.gcount();

                m_stream.next_in = m_input_buffer;
                m_stream.avail_in = data_read;
            }

            // copy the remaining m_output_leftover at the start of the stream
            memmove(m_output_buffer, m_output_buffer + m_output_sz, m_output_leftover_sz);
            uint64_t output_buffer_sz = m_buffer_capacity - m_output_leftover_sz;

            // decompress the input
            m_stream.next_out = m_output_buffer + m_output_leftover_sz;
            m_stream.avail_out = output_buffer_sz;
            int rc = inflate(&m_stream, Z_NO_FLUSH);
            if(rc != Z_OK && rc != Z_STREAM_END) ERROR("Cannot decompress the input stream: " << m_stream.msg << " (rc: " << rc << ")");
            uint64_t data_processed = m_output_leftover_sz + (output_buffer_sz - m_stream.avail_out);

            m_output_pos = 0;
            m_output_sz = (/* truncate */ data_processed / sizeof(T)) * sizeof(T);
            m_output_leftover_sz = data_processed % sizeof(T);
            assert((rc == Z_OK) || (rc == Z_STREAM_END || m_output_leftover_sz == 0));
        }

        assert(m_output_sz > 0 && "Empty output");
        assert(m_output_pos < m_output_sz && "Current output buffer depleted");
        assign(item, reinterpret_cast<T*>(m_output_buffer + m_output_pos));
        m_output_pos += sizeof(T);
        return true;
    }
};

// zlib format, non weighted
class GraphalyticsZlibReader2 : public GraphalyticsReaderBaseImpl {
private:
    GraphalyticsZlibDecompressInput<uint64_t> m_vertices; // reader for the vertex file
    GraphalyticsZlibDecompressInput<uint64_t[2]> m_edges; // reader for the edge file

public:
    GraphalyticsZlibReader2(const string& path_vertex_file, const string& path_edge_file) : m_vertices(path_vertex_file), m_edges(path_edge_file) {
        COUT_DEBUG("zlib, non weighted");
    }

    // read one edge at the time from the edge file
    bool read_edge(uint64_t& source, uint64_t& destination, double& weight) override {
        uint64_t edge[2];
        bool result = m_edges.read(&edge);
        if(!result) return false;
        source = edge[0];
        destination = edge[1];
        weight = 0;
        return true;
    }

    // read one vertex at the time from the vertex file
    bool read_vertex(uint64_t& out_vertex) override {
        return m_vertices.read(&out_vertex);
    }
};

// zlib format, weighted
class GraphalyticsZlibReader3 : public GraphalyticsReaderBaseImpl {
private:
    GraphalyticsZlibDecompressInput<uint64_t> m_vertices; // reader for the vertex file
    GraphalyticsZlibDecompressInput<uint64_t[3]> m_edges; // reader for the edge file

public:
    GraphalyticsZlibReader3(const string& path_vertex_file, const string& path_edge_file) : m_vertices(path_vertex_file), m_edges(path_edge_file) {
        COUT_DEBUG("zlib, weighted");
    }

    // read one edge at the time from the edge file
    bool read_edge(uint64_t& source, uint64_t& destination, double& weight) override {
        uint64_t edge[3];
        bool result = m_edges.read(&edge);
        if(!result) return false;
        source = edge[0];
        destination = edge[1];
        weight = reinterpret_cast<double*>(&edge)[2];
        return true;
    }

    // read one vertex at the time from the vertex file
    bool read_vertex(uint64_t& out_vertex) override {
        return m_vertices.read(&out_vertex);
    }
};


} // namespace details

/*****************************************************************************
 *                                                                           *
 *  Interface                                                                *
 *                                                                           *
 *****************************************************************************/
GraphalyticsReader::GraphalyticsReader(const std::string& path_properties) : m_random_generator(configuration().seed() + 12908478) {
    if(!common::filesystem::file_exists(path_properties)) ERROR("The given file does not exist: " << path_properties);
    string abs_path_properties = common::filesystem::absolute_path(path_properties);
    m_properties.insert({string("property-file"), abs_path_properties});
    string basedir = common::filesystem::directory(abs_path_properties);

    regex pattern { "^\\s*graph\\.([A-Za-z0-9_.-]+?)\\s*=\\s*([^#\n]+?)\\s*" };

    string name;

    COUT_DEBUG("Parsing the property file: " << path_properties);
    fstream handle { path_properties.c_str() };
    while(handle.good()){
        string line;
        getline(handle, line);

        smatch matches;
        if ( regex_match(line, matches, pattern) ) { // side effect => populate matches
            string key = matches[1];
            string value = matches[2];

            if(name.empty()){
                size_t pos = key.find('.');
                if(pos == string::npos) ERROR("Cannot parse the name of the graph (expected to find graph.<GRAPH_NAME>.<PROPERTY_NAME>), key=graph." << key);
                name = key.substr(0, pos);
                COUT_DEBUG("name=" << name);
                m_properties.insert({string("name"), name});
            }

            if( key.substr(0, name.length()) != name ){
                LOG("[GraphalyticsReader] Warning, line skipped `" << line << "', graph name does not match what expected: " << name);
                continue;
            }
            key = key.substr(name.length() +1);

            if(key == "vertex-file" || key == "edge-file") {
                if(value.empty()) ERROR("Empty path specified for the key " << key << ", line: " << line);
                if(value[0] != '/'){
                    value = common::filesystem::absolute_path(basedir + "/" + value);
                }
            } else if (key == "directed"){
                regex pattern_yes { "^\\s*(yes|true|1)\\s*$", regex_constants::icase };
                regex pattern_no { "^\\s*(no|false|0)\\s*$", regex_constants::icase };

                if(regex_search(value, pattern_yes)){
                    COUT_DEBUG("The graph is directed");
                    m_directed = true;
                } else if (regex_search(value, pattern_no)){
                    COUT_DEBUG("The graph is undirected");
                    m_directed = false;
                } else {
                    ERROR("Cannot determine whether the graph is directed or not. The property value is `" << value << "'");
                }
            } else if ( key == "compression" ){
                if(value == "zlib"){
                    m_is_compressed = true;
                } else {
                    ERROR("Compression method not supported: " << value);
                }
            } else if ( key == "edge-properties.names") {
                m_is_weighted = true;
            }

            COUT_DEBUG("key: " << key << ", value: " << value);
            m_properties.insert({key, value});
        }
    }
    handle.close();

    // check that the key vertex-file, edge-file and directed are present
    if(m_properties.find("vertex-file") == m_properties.end()) ERROR("The property `vertex-file' is not set in the property file");
    if(m_properties.find("edge-file") == m_properties.end()) ERROR("The property `edge-file' is not set in the property file");
    COUT_DEBUG("vertex-file: " << get_path_vertex_list() << ", edge-file: " << get_path_edge_list() << ", is_directed: " << is_directed() << ", is_weighted: " << is_weighted() << ", is_compressed: " << is_compressed());

    // init the reader impl. (plain or compressed)
    reset();
}

GraphalyticsReader::~GraphalyticsReader(){
    delete m_impl; m_impl = nullptr;
}

void GraphalyticsReader::reset(){
    delete m_impl; m_impl = nullptr;

    if(!is_compressed()){
        m_impl = new GraphalyticsPlainReader{ get_path_vertex_list(), get_path_edge_list(), is_weighted() };
    } else {
        if(!is_weighted()){
            m_impl = new GraphalyticsZlibReader2( get_path_vertex_list(), get_path_edge_list() );
        } else {
            m_impl = new GraphalyticsZlibReader3( get_path_vertex_list(), get_path_edge_list() );
        }
    }
}

string GraphalyticsReader::get_property(const std::string& key) const {
    auto item = m_properties.find(key);
    return (item == m_properties.end()) ? "" : item->second;
}

string GraphalyticsReader::get_path_vertex_list() const {
    auto item = m_properties.find("vertex-file");
    if(item == m_properties.end()) ERROR("The property `vertex-file' is not set");
    return item->second;
}

string GraphalyticsReader::get_path_edge_list() const {
    auto item = m_properties.find("edge-file");
    if(item == m_properties.end()) ERROR("The property `edge-file' is not set");
    return item->second;
}

bool GraphalyticsReader::is_directed() const {
    return m_directed;
}

bool GraphalyticsReader::is_weighted() const {
    return m_is_weighted;
}

bool GraphalyticsReader::is_compressed() const {
    return m_is_compressed;
}

void GraphalyticsReader::set_emit_directed_edges(bool value){
    m_emit_directed_edges = value;
}

bool GraphalyticsReader::read(graph::WeightedEdge& edge){
    return read_edge(edge);
}

bool GraphalyticsReader::read_edge(graph::WeightedEdge& edge){
    if(!is_directed() && m_emit_directed_edges && !m_last_reported){
        std::swap(m_last_source, m_last_destination);
        m_last_reported = true;
    } else {
        bool has_result = m_impl->read_edge(m_last_source, m_last_destination, m_last_weight);
        if(!has_result) return false;

        if(!is_weighted()){
            uniform_real_distribution<double> distribution{0, configuration().max_weight()}; // in [0, max_weight)
            m_last_weight = distribution(m_random_generator);
            if(m_last_weight == 0.0) m_last_weight = configuration().max_weight(); // in (0, max_weight]
        }

        m_last_reported = false;
    }

    edge.m_source = m_last_source;
    edge.m_destination = m_last_destination;
    edge.m_weight = m_last_weight;
    COUT_DEBUG("edge parsed: " << edge);

    return true;
}

bool GraphalyticsReader::read_vertex(uint64_t& out_vertex){
    return m_impl->read_vertex(out_vertex);
}


} // namespace

