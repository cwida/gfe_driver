#include "gtest/gtest.h"

#include "lib/common/static_index.hpp"

#include <iostream>
#include <memory>
#include <unistd.h> // sleep
#include <utility>

using namespace std;
using namespace common;

TEST(StaticIndex, only_root){
    StaticIndex<int64_t> index(/* node size */ 4, /* number of keys */ 3);
    index.set_separator_key(0, 10);
    index.set_separator_key(1, 20);
    index.set_separator_key(2, 30);

    ASSERT_EQ(index.get_separator_key(0), 10);
    ASSERT_EQ(index.get_separator_key(1), 20);
    ASSERT_EQ(index.get_separator_key(2), 30);

    ASSERT_EQ(index.find_lte(5), 0);
    ASSERT_EQ(index.find_lte(10), 0);
    ASSERT_EQ(index.find_lte(15), 0);
    ASSERT_EQ(index.find_lte(20), 1);
    ASSERT_EQ(index.find_lte(25), 1);
    ASSERT_EQ(index.find_lte(30), 2);
    ASSERT_EQ(index.find_lte(35), 2);

    ASSERT_EQ(index.find_lt(5), 0);
    ASSERT_EQ(index.find_lt(10), 0);
    ASSERT_EQ(index.find_lt(15), 0);
    ASSERT_EQ(index.find_lt(20), 0);
    ASSERT_EQ(index.find_lt(25), 1);
    ASSERT_EQ(index.find_lt(30), 1);
    ASSERT_EQ(index.find_lt(35), 2);

    ASSERT_EQ(index.find_lte_last(5), 0);
    ASSERT_EQ(index.find_lte_last(10), 0);
    ASSERT_EQ(index.find_lte_last(15), 0);
    ASSERT_EQ(index.find_lte_last(20), 1);
    ASSERT_EQ(index.find_lte_last(25), 1);
    ASSERT_EQ(index.find_lte_last(30), 2);
    ASSERT_EQ(index.find_lte_last(35), 2);
}

TEST(StaticIndex, height2){
    StaticIndex<int> index(/* node size */ 4, /* number of keys */ 7);
    for(int i = 0; i < 7; i++){ index.set_separator_key(i, (i+1) * 10); } // 10, 20, 30, 40, 50, 60, 70

    // check
    for(int i = 0; i < 7; i++) { ASSERT_EQ(index.get_separator_key(i), (i+1) *10); }

    ASSERT_EQ(index.find_lte(5), 0);
    for(int64_t key = 10; key <= 75; key+=5){
        auto segment_id = index.find_lte(key);
        ASSERT_EQ(segment_id, (key / 10) -1);
    }

    for(int64_t key = 5; key <= 20; key += 5){ ASSERT_EQ(index.find_lt(key), 0); }
    for(int64_t key = 25; key <= 75; key += 5){
        int64_t expected_segment = ((key -1) / 10) -1;
        ASSERT_EQ(index.find_lt(key), expected_segment);
    }
    for(int64_t key = 75; key >= 5; key -= 5){
        int64_t expected_segment = max<int64_t>(( key - 10 ) / 10, 0);
        ASSERT_EQ(index.find_lte_last(key), expected_segment);
    }
}

TEST(StaticIndex, full_tree){
    constexpr size_t num_keys = 64;

    StaticIndex<int> index(/* node size */ 4, /* number of keys */ num_keys);
    for(int i = 0; i < num_keys; i++){
        index.set_separator_key(i, (i+1) * 10);
        for(int j = 0; j <= i; j++) {
            ASSERT_EQ(index.get_separator_key(j), (j+1) * 10);
        }
    } // 10, 20, 30, 40, 50, 60, 70, etc.

    // check
    for(int i = 0; i < num_keys; i++){
        ASSERT_EQ(index.find_lte((i+1) * 10 -1), max(i -1, 0));
        ASSERT_EQ(index.find_lte((i+1) * 10), i);
        ASSERT_EQ(index.find_lte((i+1) * 10 +1), i);
    }
}

TEST(StaticIndex, height5){
    constexpr size_t num_keys = 366;

    StaticIndex<int> index(/* node size */ 4, /* number of keys */ num_keys);
    for(int i = 0; i < num_keys; i++){
        index.set_separator_key(i, (i+1) * 10);
        for(int j = 0; j <= i; j++) {
            ASSERT_EQ(index.get_separator_key(j), (j+1) * 10);
        }
    } // 10, 20, 30, 40, 50, 60, 70, etc.

    // check
    for(int i = 0; i < num_keys; i++){
        ASSERT_EQ(index.find_lte((i+1) * 10 -1), max(i -1, 0));
        ASSERT_EQ(index.find_lte((i+1) * 10), i);
        ASSERT_EQ(index.find_lte((i+1) * 10 +1), i);
    }
}

TEST(StaticIndex, height6){
    constexpr size_t num_keys = 4000;

    StaticIndex<int> index(/* node size */ 5, /* number of keys */ num_keys);
    ASSERT_EQ(index.height(), 6);
    for(int i = 0; i < num_keys; i++){
        index.set_separator_key(i, (i+1) * 10);
    } // 10, 20, 30, 40, 50, 60, 70, etc.

    for(int i = 0; i < num_keys; i++) {
        ASSERT_EQ(index.get_separator_key(i), (i+1) * 10);
    }

    // check
    for(int i = 0; i < num_keys; i++){
        ASSERT_EQ(index.find_lte((i+1) * 10 -1), max(i -1, 0));
        ASSERT_EQ(index.find_lte((i+1) * 10), i);
        ASSERT_EQ(index.find_lte((i+1) * 10 +1), i);
    }
}

TEST(StaticIndex, duplicates){
    StaticIndex<int> index(/* node size */ 4, /* number of keys */ 7);
    index.set_separator_key(0, 5);
    index.set_separator_key(1, 5);
    index.set_separator_key(2, 10);
    index.set_separator_key(3, 10);
    index.set_separator_key(4, 15);
    index.set_separator_key(5, 20);
    index.set_separator_key(6, 25);

    ASSERT_EQ(index.find_lt(4), 0);
    ASSERT_EQ(index.find_lt(5), 0);
    ASSERT_EQ(index.find_lt(9), 1);
    ASSERT_EQ(index.find_lt(10), 1);
    ASSERT_EQ(index.find_lt(11), 3);

    ASSERT_EQ(index.find_lte_first(4), 0);
    ASSERT_EQ(index.find_lte_first(5), 0);
    ASSERT_EQ(index.find_lte_first(6), 1);
    ASSERT_EQ(index.find_lte_first(9), 1);
    ASSERT_EQ(index.find_lte_first(10), 2);
    ASSERT_EQ(index.find_lte_first(11), 3);
    ASSERT_EQ(index.find_lte_first(14), 3);
    ASSERT_EQ(index.find_lte_first(15), 4);
    ASSERT_EQ(index.find_lte_first(24), 5);
    ASSERT_EQ(index.find_lte_first(25), 6);
    ASSERT_EQ(index.find_lte_first(26), 6);

    ASSERT_EQ(index.find_lte_last(4), 0);
    ASSERT_EQ(index.find_lte_last(5), 1);
    ASSERT_EQ(index.find_lte_last(6), 1);
    ASSERT_EQ(index.find_lte_last(9), 1);
    ASSERT_EQ(index.find_lte_last(10), 3);
    ASSERT_EQ(index.find_lte_last(11), 3);
    ASSERT_EQ(index.find_lte_last(14), 3);
    ASSERT_EQ(index.find_lte_last(15), 4);
    ASSERT_EQ(index.find_lte_last(24), 5);
    ASSERT_EQ(index.find_lte_last(25), 6);
    ASSERT_EQ(index.find_lte_last(26), 6);
}
