#include "gtest/gtest.h"

#include <algorithm>
#include <cinttypes>
#include "lib/common/sorting.hpp"

using namespace std;
using namespace common;

TEST(Sorting, sanity){
    uint64_t array[] = {50, 60, 30, 10, 40, 80, 90, 70, 20};
    const uint64_t array_sz = sizeof(array) / sizeof(array[0]);

    details::sorting::implementation(array, array_sz, std::less<uint64_t>(), /* num samples */ 2);

    for(uint64_t i = 0; i < array_sz; i++){
        ASSERT_EQ(array[i], (i +1) * 10);
    }
}

TEST(Sorting, empty){
    uint64_t array[] = {};
    const uint64_t array_sz = sizeof(array) / sizeof(array[0]);

    details::sorting::implementation(array, array_sz, std::less<uint64_t>(), /* num samples */ 2);
}

TEST(Sorting, single_element){
    uint64_t array[] = {10};
    const uint64_t array_sz = sizeof(array) / sizeof(array[0]);

    details::sorting::implementation(array, array_sz, std::less<uint64_t>(), /* num samples */ 2);

    ASSERT_EQ(array[0], 10);
}

TEST(Sorting, two_elements){
    uint64_t array[] = {20, 10};
    const uint64_t array_sz = sizeof(array) / sizeof(array[0]);

    details::sorting::implementation(array, array_sz, std::less<uint64_t>(), /* num samples */ 2);

    ASSERT_EQ(array[0], 10);
    ASSERT_EQ(array[1], 20);
}

TEST(Sorting, reverse){
    uint64_t array[] = {50, 60, 30, 10, 40, 80, 90, 70, 20};
    const uint64_t array_sz = sizeof(array) / sizeof(array[0]);

    details::sorting::implementation(array, array_sz, std::greater<uint64_t>(), /* num samples */ 2);

    for(uint64_t i = 0; i < array_sz; i++){
        ASSERT_EQ(array[array_sz - 1 - i], (i +1) * 10);
    }
}


