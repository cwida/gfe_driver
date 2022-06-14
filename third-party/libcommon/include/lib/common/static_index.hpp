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

#ifndef COMMON_STATIC_INDEX_HPP
#define COMMON_STATIC_INDEX_HPP

#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>

namespace common {

/**
 * A static index is an index a fixed number of entries. The key can be any fixed-length and trivially
 * copyable data table (int, double), while the values are implicitly the integers in [0, number of entries).
 * To change the number of entries,  the whole index needs to be rebuilt (#rebuild(N)) providing the
 * new size of the size.
 *
 * The node size B is determined on initialisation. A node size B actually requires B -1 slots
 * in terms of space, so it is recommended to set B to a power of 2 + 1 (e.g. 65) to fully
 * exploit aligned accesses to the cache.
 *
 * This class is not thread-safe.
 */
template<typename KeyType>
class StaticIndex {
    const uint16_t m_node_size; // number of keys per node
    int16_t m_height; // the height of this tree
    int32_t m_capacity; // the number of segments/keys in the tree
    KeyType* m_keys; // the container of the keys
    KeyType m_key_minimum; // the minimum stored in the tree

    /**
     * Keep track of the cardinality and the height of the rightmost subtrees
     */
    struct RightmostSubtreeInfo {
        uint16_t m_root_sz; // the number of elements in the root
        uint16_t m_right_height; // the height of the rightmost subtree
    };
    constexpr static uint64_t m_rightmost_sz = 8;
    RightmostSubtreeInfo m_rightmost[m_rightmost_sz];

protected:
    // Retrieve the slot associated to the given entry
    KeyType* get_slot(uint64_t position) const;

    // Dump the content of the given subtree
    void dump_subtree(std::ostream& out, KeyType* root, int height, bool rightmost, KeyType fence_min, KeyType fence_max, bool* integrity_check) const;
    static void dump_tabs(std::ostream& out, size_t depth);

    // Generic implementation of the method #find
    template <typename Fn>
    uint64_t traverse_tree(KeyType key, Fn fn) const noexcept;

public:
    /**
     * Initialise the AB-Tree with the given node size and capacity
     */
    StaticIndex(uint64_t node_size = 65, uint64_t num_entries = 1);

    /**
     * Destructor
     */
    ~StaticIndex();

    /**
     * Rebuild the tree to contain `num_entries'
     */
    void rebuild(uint64_t num_entries);

    /**
     * Set the separator key associated to the given entry
     */
    void set_separator_key(uint64_t position, KeyType key);

    /**
     * Get the separator key associated to the given entry.
     * Used only for the debugging purposes.
     */
    KeyType get_separator_key(uint64_t position) const;

    /**
     * Return the entry related to the largest key smaller than the given key
     */
    uint64_t find_lt(KeyType key) const noexcept;

    /**
     * Return an entry that contains the given key. If there are no repetitions in the indexed data structure,
     * this will be the only candidate entry for the given key.
     */
    uint64_t find_lte(KeyType key) const noexcept;

    /**
     * Return the first entry such that its key is smaller or equal than the given key.
     */
    uint64_t find_lte_first(KeyType key) const noexcept;

    /**
     * Return the last entry such that its key is smaller or equal than key. This is useful
     * when the index is referring to the leaves of a B-Tree/PMA, then the method returns
     * the last leaf/segment that may contain the end key of a range search [start, end].
     */
    uint64_t find_lte_last(KeyType key) const noexcept;

    /**
     * Retrieve the minimum key stored in the tree
     */
    KeyType minimum() const noexcept;

    /**
     * Retrieve the height of the current static tree
     */
    int height() const noexcept;

    /**
     * Retrieve the block size of each node in the tree
     */
    int64_t node_size() const noexcept;

    /**
     * Retrieve the memory footprint of this index, in bytes
     */
    size_t memory_footprint() const;

    /**
     * Dump the fields of the index
     */
    void dump(std::ostream& out, bool* integrity_check = nullptr) const;
    void dump() const;
};

/**
 * Dump the the content of the static index
 */
template<typename KeyType>
std::ostream& operator<<(std::ostream& out, const StaticIndex<KeyType>& index);


/*****************************************************************************
 *                                                                           *
 *   Implementation details                                                  *
 *                                                                           *
 *****************************************************************************/

template<typename KeyType>
StaticIndex<KeyType>::StaticIndex(uint64_t node_size, uint64_t num_entries) :
        m_node_size(node_size), m_height(0), m_capacity(0), m_keys(nullptr), m_key_minimum(std::numeric_limits<KeyType>::max()) {
    if(node_size > (uint64_t) std::numeric_limits<uint16_t>::max()){ throw std::invalid_argument("Invalid node size: too big"); }
    rebuild(num_entries);
}

template<typename KeyType>
StaticIndex<KeyType>::~StaticIndex(){
    free(m_keys); m_keys = nullptr;
}

template<typename KeyType>
int64_t StaticIndex<KeyType>::node_size() const noexcept {
    // cast to int64_t
    return m_node_size;
}

template<typename KeyType>
void StaticIndex<KeyType>::rebuild(uint64_t N){
    if(N == 0) throw std::invalid_argument("Invalid number of keys: 0");
    int height = ceil( log2(N) / log2(node_size()) );
    if(height > static_cast<int>(m_rightmost_sz)){ throw std::invalid_argument("Invalid number of keys/segments: too big"); }
    uint64_t tree_sz = pow(node_size(), height) -1; // don't store the minimum, segment 0

    if(height != m_height){
        free(m_keys); m_keys = nullptr;
        int rc = posix_memalign((void**) &m_keys, /* alignment */ 64,  /* size */ tree_sz * sizeof(KeyType));
        if(rc != 0) { throw std::bad_alloc(); }
        m_height = height;
    }
    m_capacity = N;

    // set the height of all rightmost subtrees
    while(height > 0){
        uint64_t subtree_sz = pow(node_size(), height -1);
        m_rightmost[height - 1].m_root_sz = (N -1) / subtree_sz;
        assert(m_rightmost[height -1].m_root_sz > 0);
        uint64_t rightmost_subtree_sz = (N -1) % subtree_sz;
        int rightmost_subtree_height = 0;
        if(rightmost_subtree_sz > 0){
            rightmost_subtree_sz += 1; // with B-1 keys we index B entries
            rightmost_subtree_height = ceil(log2(rightmost_subtree_sz) / log2(m_node_size));
        }
        m_rightmost[height -1].m_right_height = rightmost_subtree_height;

        // next subtree
        N = rightmost_subtree_sz;
        height = rightmost_subtree_height;
    }
}

template<typename KeyType>
int StaticIndex<KeyType>::height() const noexcept {
    return m_height;
}

template<typename KeyType>
size_t StaticIndex<KeyType>::memory_footprint() const {
    return (pow(node_size(), height()) -1) * sizeof(KeyType);
}

template<typename KeyType>
KeyType* StaticIndex<KeyType>::get_slot(uint64_t entry_id) const {
    assert(entry_id > 0 && "The segment 0 is not explicitly stored");
    assert(entry_id < static_cast<uint64_t>(m_capacity) && "Invalid slot");

    KeyType* __restrict base = m_keys;
    int64_t offset = entry_id;
    int height = m_height;
    bool rightmost = true; // this is the rightmost subtree
    int64_t subtree_sz = pow(node_size(), height -1);

    while(height > 0){
        int64_t subtree_id = offset / subtree_sz;
        int64_t modulo = offset % subtree_sz;

        if(modulo == 0){ // found, this is an internal node
            assert(subtree_id > 0 && "Otherwise this would have been an internal element on an ancestor");
            return base + subtree_id -1;
        }

        // traverse the children
        base += (node_size() -1) + subtree_id * (subtree_sz -1);
        offset -= subtree_id * subtree_sz;

        // is this the rightmost subtree ?
        rightmost = rightmost && (subtree_id >= m_rightmost[height -1].m_root_sz);
        if(rightmost){
            height = m_rightmost[height -1].m_right_height;
            subtree_sz = pow(node_size(), height -1);
        } else {
            height --;
            subtree_sz /= node_size();
        }
    }

    return base + offset;
}

template<typename KeyType>
void StaticIndex<KeyType>::set_separator_key(uint64_t position, KeyType key){
    if(position == 0) {
        m_key_minimum = key;
    } else {
        get_slot(position)[0] = key;
    }

    assert(get_separator_key(position) == key);
}

template<typename KeyType>
KeyType StaticIndex<KeyType>::get_separator_key(uint64_t position) const {
    if(position == 0)
        return m_key_minimum;
    else
        return get_slot(position)[0];
}

template <typename KeyType>
template <typename Fn>
uint64_t StaticIndex<KeyType>::traverse_tree(KeyType key, Fn fn) const noexcept{
    KeyType* __restrict base = m_keys;
    int64_t offset = 0;
    int height = m_height;
    bool rightmost = true; // this is the rightmost subtree
    int64_t subtree_sz = pow(node_size(), height -1);

    while(height > 0){
        uint64_t root_sz = (rightmost) ? m_rightmost[height -1].m_root_sz : node_size() -1; // full
        uint64_t subtree_id = fn(base, root_sz, key);

        base += (node_size() -1) + subtree_id * (subtree_sz -1);
        offset += subtree_id * subtree_sz;

        // similar to #get_slot
        rightmost = rightmost && (subtree_id >= m_rightmost[height -1].m_root_sz);
        if(rightmost){
            height = m_rightmost[height -1].m_right_height;
            subtree_sz = pow(node_size(), height -1);
        } else {
            height --;
            subtree_sz /= node_size();
        }
    }

    return offset;
}

template<typename KeyType>
uint64_t StaticIndex<KeyType>::find_lt(KeyType key) const noexcept {
    if(key < m_key_minimum) return 0; // easy!

    return traverse_tree(key, [](KeyType* __restrict node, uint64_t node_sz, KeyType key){
        uint64_t position = 0;
        while(position < node_sz && node[position] < key) position++;
        return position;
    });
}

template<typename KeyType>
uint64_t StaticIndex<KeyType>::find_lte(KeyType key) const noexcept {
    if(key <= m_key_minimum) return 0; // easy!

    return traverse_tree(key, [](KeyType* __restrict node, uint64_t node_sz, KeyType key){
        uint64_t position = 0;
        while(position < node_sz && node[position] <= key) position++;
        return position;
    });
}

template<typename KeyType>
uint64_t StaticIndex<KeyType>::find_lte_first(KeyType key) const noexcept {
    if(key <= m_key_minimum) return 0; // easy!

    return traverse_tree(key, [](KeyType* __restrict node, uint64_t node_sz, KeyType key){
        uint64_t position = 0;
        while(position < node_sz && node[position] < key) { position++; }
        if(position < node_sz && node[position] == key) { position++; }
        return position;
    });
}

template<typename KeyType>
uint64_t StaticIndex<KeyType>::find_lte_last(KeyType key) const noexcept {
    if(key < m_key_minimum) return 0; // easy!

    return traverse_tree(key, [](KeyType* __restrict node, uint64_t node_sz, KeyType key){
        uint64_t position = node_sz;
        while(position > 0 && key < node[position -1]) position--;
        return position;
    });
}

template<typename KeyType>
KeyType StaticIndex<KeyType>::minimum() const noexcept {
    return m_key_minimum;
}

template<typename KeyType>
void StaticIndex<KeyType>::dump_tabs(std::ostream& out, size_t depth){
    using namespace std;

    auto flags = out.flags();
    out << std::setw((depth-1) * 2 + 5) << setfill(' ') << ' ';
    out.setf(flags);
}

template<typename KeyType>
void StaticIndex<KeyType>::dump_subtree(std::ostream& out, KeyType* root, int height, bool rightmost, KeyType fence_min, KeyType fence_max, bool* integrity_check) const {
    using namespace std;
    if(height <= 0) return; // base case

    int depth = m_height - height +1;
    int64_t root_sz = (rightmost) ? m_rightmost[height -1].m_root_sz : node_size() -1; // full
    int64_t subtree_sz = pow(node_size(), height -1);

    // preamble
    auto flags = out.flags();
    if(depth > 1) out << ' ';
    out << setw((depth -1) * 2) << setfill(' '); // initial padding
    out << "[" << setw(2) << setfill('0') << depth << "] ";
    out << "offset: " << root - m_keys << ", root size: " << root_sz << ", fence keys (interval): [" << fence_min << ", " << fence_max << "]\n";
    out.setf(flags);

    dump_tabs(out, depth);
    out << "keys: ";
    for(size_t i = 0; i < root_sz; i++){
        if(i > 0) out << ", ";
        out << (i+1) << " => k:" << (i+1) * subtree_sz << ", v:" << root[i];

        if(i == 0 && root[0] < fence_min){
            out << " (ERROR: smaller than the min fence key: " << fence_min << ")";
            if(integrity_check) *integrity_check = false;
        }
        if(i > 0 && root[i] < root[i-1]){
            out << " (ERROR: sorted order not respected: " << root[i-1] << " > " << root[i] << ")";
            if(integrity_check) *integrity_check = false;
        }
        if(i == root_sz -1 && root[i] > fence_max){
            out << " (ERROR: greater than the max fence key: " << fence_max << ")";
            if(integrity_check) *integrity_check = false;
        }

    }
    out << "\n";

    if(height > 1) { // internal node?
        KeyType* base = root + node_size() -1;

        dump_tabs(out, depth);
        out << "offsets: ";
        for(size_t i = 0; i <= root_sz; i++){
          if(i > 0) out << ", ";
          out << i << ": " << (base + i * (subtree_sz -1)) - m_keys;
        }
        out << "\n";

        // recursively dump the children
        for(size_t i = 0; i < root_sz; i++){
            KeyType fmin = (i == 0) ? fence_min : root[i-1];
            KeyType fmax = root[i];

            dump_subtree(out, base + (i* (subtree_sz -1)), height -1, false, fmin, fmax, integrity_check);
        }

        // dump the rightmost subtree
        dump_subtree(out, base + root_sz * (subtree_sz -1), m_rightmost[height -1].m_right_height, rightmost, root[root_sz -1], fence_max, integrity_check);
    }
}

template<typename KeyType>
void StaticIndex<KeyType>::dump(std::ostream& out, bool* integrity_check) const {
    out << "[Index] block size: " << node_size() << ", height: " << height() <<
            ", capacity (number of entries indexed): " << m_capacity << ", minimum: " << minimum() << "\n";

    if(m_capacity > 1)
        dump_subtree(out, m_keys, height(), true, m_key_minimum, std::numeric_limits<KeyType>::max(), integrity_check);
}

template<typename KeyType>
void StaticIndex<KeyType>::dump() const {
    dump(std::cout);
}

template<typename KeyType>
std::ostream& operator<<(std::ostream& out, const StaticIndex<KeyType>& index){
    index.dump(out);
    return out;
}

} // namespace common

#endif /* COMMON_STATIC_INDEX_HPP */
