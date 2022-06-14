#include "gtest/gtest.h"

#include <cinttypes>
#include <climits>
#include <cstdio>
#include <random>
#include "lib/common/profiler.hpp"

using namespace std;
using namespace common;

volatile uint64_t cache_faults_count = 0;

TEST(Profiler, cache_faults){
    constexpr int64_t A_sz = (1ull << 20); // 1M
    int64_t* A = new int64_t[A_sz];

    // init
    mt19937_64 random{1ull};
    CachesProfiler c;
    c.start();
    A[0] = 0;
    for(int64_t i = 1; i < A_sz; i++){
        int64_t j = uniform_int_distribution<int64_t>{0, i-1}(random);
        A[i] = A[j];
        A[j] = i; // ptr to A[i];
    }
    c.stop();

    cout << "Init, sequential accesses, faults: " << c.snapshot() << endl;

    // run
    c.start();
    int64_t index = 0;
    for(int i = 1; i < (1ull << 20); i++){
        index = A[index];
        cache_faults_count += index;
    }
    c.stop();
    cout << "Run, random accesses, faults: " << c.snapshot() << endl;

    delete[] A; A = nullptr;
}

TEST(Profiler, software_events){
    cout << "\n";

    // Create a temporary file
    char path[PATH_MAX];
    strcpy(path, "/tmp/test_profiler_XXXXXX");
    int fd = mkstemp(path);
    ASSERT_NE(fd, -1); // cannot create a temporary file
    cout << "Temporary file: " << path << endl;
    FILE* file = fdopen(fd, "w");
    uint64_t cardinality = 1024 * 1024; // 8 MB
    for(uint64_t i = 0; i < cardinality; i++){
        size_t rc = fwrite(&i, sizeof(i), 1, file);
        ASSERT_EQ(rc, 1); // it should have written exactly one element
    }
    int rc = fclose(file); // it also closes the underlying file descriptor
    ASSERT_EQ(rc, 0);

    SoftwareEventsProfiler profiler;
    profiler.start();

    // Read the content of the temporary file just created
    file = fopen(path, "r");
    ASSERT_NE(file, nullptr);
    for(uint64_t i = 0; i < cardinality; i++){
        uint64_t value = 0;
        size_t count = fread(&value, sizeof(value), 1, file);
        ASSERT_EQ(count, 1);
        ASSERT_EQ(value, i);
    }

    rc = fclose(file);
    ASSERT_EQ(rc, 0);

    profiler.stop();
    cout << "Software events: " << profiler.snapshot() << endl;
}