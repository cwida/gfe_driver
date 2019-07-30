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

ExecStatistics::ExecStatistics(std::vector<uint64_t>& trials) : m_num_trials(trials.size()){
    if(m_num_trials == 0){
        memset((void*) m_percentile, 0, sizeof(m_percentile));
    } else {
        // compute mean/std.dev/min/max
        uint64_t sum = 0;
        uint64_t sum2 = 0;
        uint64_t vmin = std::numeric_limits<uint64_t>::max();
        uint64_t vmax = 0;
        for(uint64_t value: trials){
            sum += value;
            sum2 += value * value;
            vmin = min(value, vmin);
            vmax = max(value, vmax);
        }

        m_mean = sum / m_num_trials;
        m_stddev = (sum2 / m_num_trials) - (m_mean * m_mean);
        m_min = vmin;
        m_max = vmax;

        // compute the percentiles
        sort(begin(trials), end(trials));
        for(uint64_t i = 1; i <= 99; i++){
            uint64_t pos = (i * m_num_trials) / 100; // i : 100 = pos : num_trials
            m_percentile[i -1] = trials[pos > 0 ? (pos -1) : 0] ;
        }

        // special case for the median
        if(m_num_trials % 2 == 0){
            m_percentile[49] = (trials[m_num_trials /2] + trials[(m_num_trials /2)-1]) /2;
        }
    }
}


void ExecStatistics::save(const std::string& name){
    assert(configuration().db() != nullptr);

    auto store = configuration().db()->add("statistics");
    store.add("type", name);
    store.add("num_trials", m_num_trials);
    store.add("mean", m_mean);
    store.add("median", m_percentile[49]);
    store.add("stddev", m_stddev);
    store.add("min", m_min);
    store.add("max", m_max);
    for(uint64_t i = 0; i < 99; i++){
        stringstream ss_field_name;
        ss_field_name << "percentile_";
        ss_field_name << (i+1);
        string field_name = ss_field_name.str();
        store.add(field_name, m_percentile[i]);
    }
}

std::ostream& operator<<(std::ostream& out, const ExecStatistics& stats){
    out << "N: " << stats.m_num_trials << ", mean: " << stats.m_mean << ", median: " << stats.m_percentile[49] << ", "
            << "std. dev.: " << stats.m_stddev << ", min: " << stats.m_min << ", max: " << stats.m_max << ", perc 90: " << stats.m_percentile[89] << ", "
            << "perc 95: " << stats.m_percentile[94] << ", perc 99: " << stats.m_percentile[98] << "]";

    return out;
}

} // namespace


