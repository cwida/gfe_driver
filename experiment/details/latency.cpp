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

#include "latency.hpp"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <limits>
#include "common/database.hpp"
#include "common/quantity.hpp"
#include "configuration.hpp"

using namespace common;
using namespace std;

namespace gfe::experiment::details {

static uint64_t get_percentile(uint64_t* __restrict A, uint64_t A_sz, uint64_t index){
    uint64_t pos = (index * A_sz) / 100; // i : 100 = pos : num_operations
    return A[pos > 0 ? (pos -1) : 0];
}

LatencyStatistics LatencyStatistics::compute_statistics(uint64_t* arr_latencies_nanosecs, uint64_t arr_latencies_sz){
    LatencyStatistics instance;
    if(arr_latencies_sz == 0) return instance;

    // compute mean/std.dev/min/max
    uint64_t sum = 0;
    uint64_t sum2 = 0;
    uint64_t vmin = numeric_limits<uint64_t>::max();
    uint64_t vmax = 0;
    for(uint64_t i = 0; i < arr_latencies_sz; i++){
        uint64_t value = arr_latencies_nanosecs[i];

        sum += value;
        sum2 += value * value;
        vmin = min(value, vmin);
        vmax = max(value, vmax);
    }

    instance.m_num_operations = arr_latencies_sz;
    instance.m_mean = sum / arr_latencies_sz;
    instance.m_stddev = (sum2 / arr_latencies_sz) - (instance.m_mean * instance.m_mean);
    instance.m_min = vmin;
    instance.m_max = vmax;

    // compute the percentiles
    sort(arr_latencies_nanosecs, arr_latencies_nanosecs + arr_latencies_sz);
    instance.m_percentile90 = get_percentile(arr_latencies_nanosecs, arr_latencies_sz, 90);
    instance.m_percentile95 = get_percentile(arr_latencies_nanosecs, arr_latencies_sz, 95);
    instance.m_percentile97 = get_percentile(arr_latencies_nanosecs, arr_latencies_sz, 97);
    instance.m_percentile99 = get_percentile(arr_latencies_nanosecs, arr_latencies_sz, 99);

    if(arr_latencies_sz % 2 == 0){
        instance.m_median = (arr_latencies_nanosecs[arr_latencies_sz /2] + arr_latencies_nanosecs[(arr_latencies_sz /2)-1]) /2;
    } else {
        instance.m_median = arr_latencies_nanosecs[arr_latencies_sz /2];
    }

    return instance;
}

chrono::nanoseconds LatencyStatistics::mean() const {
    return chrono::nanoseconds(m_mean);
}

chrono::nanoseconds LatencyStatistics::percentile90() const {
    return chrono::nanoseconds(m_percentile90);
}

chrono::nanoseconds LatencyStatistics::percentile99() const {
    return chrono::nanoseconds(m_percentile99);
}

void LatencyStatistics::save(const std::string& name){
    assert(configuration().db() != nullptr);

    auto store = configuration().db()->add("latencies");
    store.add("type", name);
    store.add("num_operations", m_num_operations);
    store.add("mean", m_mean);
    store.add("median", m_median);
    store.add("stddev", m_stddev);
    store.add("min", m_min);
    store.add("max", m_max);
    store.add("p90", m_percentile90);
    store.add("p95", m_percentile95);
    store.add("p97", m_percentile97);
    store.add("p99", m_percentile99);
}

static DurationQuantity _D(uint64_t value){
    return DurationQuantity{chrono::nanoseconds(value)};
}

std::ostream& operator<<(std::ostream& out, const LatencyStatistics& stats){
    out << "N: " << stats.m_num_operations << ", mean: " << _D(stats.m_mean) << ", median: " << _D(stats.m_median) << ", "
            << "std. dev.: " << _D(stats.m_stddev) << ", min: " << _D(stats.m_min) << ", max: " << _D(stats.m_max) << ", "
            << "perc 90: " << _D(stats.m_percentile90) << ", perc 95: " << _D(stats.m_percentile95) << ", "
            << "perc 97: " << _D(stats.m_percentile97) << ", perc 99: " << _D(stats.m_percentile99);
    return out;
}

} // namespace
