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

#include "statistics.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <limits>

#include "common/database.hpp"
#include "common/error.hpp"
#include "configuration.hpp"

using namespace std;

namespace experiment {

/*****************************************************************************
 *                                                                           *
 *  Debug                                                                    *
 *                                                                           *
 *****************************************************************************/
#define DEBUG
#define COUT_DEBUG_FORCE(msg) { std::cout << "[Statistics @ " << __FUNCTION__ << "] " << msg << std::endl; }
#if defined(DEBUG)
    #define COUT_DEBUG(msg) COUT_DEBUG_FORCE(msg)
#else
    #define COUT_DEBUG(msg)
#endif



/*****************************************************************************
 *                                                                           *
 *  ExecStatistics                                                           *
 *                                                                           *
 ****************************************************************************/

ExecStatistics::ExecStatistics(const std::vector<int64_t>& trials) : m_num_trials(trials.size()){
    vector<uint64_t> values;
    values.reserve(trials.size());

    // compute mean/std.dev/min/max
    uint64_t sum = 0;
    uint64_t sum2 = 0;
    uint64_t vmin = std::numeric_limits<uint64_t>::max();
    uint64_t vmax = 0;
    for(int64_t t: trials){
        if(t < 0){ m_num_timeouts++; continue; }
        uint64_t value = static_cast<uint64_t>(t);
        values.push_back(value);

        sum += value;
        sum2 += value * value;
        vmin = min(value, vmin);
        vmax = max(value, vmax);
    }

    const uint64_t num_values = values.size();
    if(num_values > 0){
        m_mean = sum / num_values;
        m_stddev = (sum2 / num_values) - (m_mean * m_mean);
        m_min = vmin;
        m_max = vmax;

        // compute the percentiles
        sort(begin(values), end(values));
        m_percentile90 = get_percentile(values, 90);
        m_percentile95 = get_percentile(values, 95);
        m_percentile97 = get_percentile(values, 97);
        m_percentile99 = get_percentile(values, 99);

        if(num_values % 2 == 0){
            m_median = (values[num_values /2] + values[(num_values /2)-1]) /2;
        }
    }

}

uint64_t ExecStatistics::get_percentile(const std::vector<uint64_t>& values_sorted, uint64_t index){
    uint64_t pos = (index * values_sorted.size()) / 100; // i : 100 = pos : num_trials
    return values_sorted[pos > 0 ? (pos -1) : 0];
}

void ExecStatistics::save(const std::string& name){
    assert(configuration().db() != nullptr);

    auto store = configuration().db()->add("statistics");
    store.add("type", name);
    store.add("num_trials", m_num_trials);
    store.add("num_timeouts", m_num_timeouts);
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

std::ostream& operator<<(std::ostream& out, const ExecStatistics& stats){
    out << "N: " << stats.m_num_trials << ", mean: " << stats.m_mean << ", median: " << stats.m_median << ", "
            << "std. dev.: " << stats.m_stddev << ", min: " << stats.m_min << ", max: " << stats.m_max << ", perc 90: " << stats.m_percentile90 << ", "
            << "perc 95: " << stats.m_percentile95 << ", perc 99: " << stats.m_percentile99 << ", num timeouts: " << stats.m_num_timeouts << "]";

    return out;
}

} // namespace


