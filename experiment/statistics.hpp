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
#include <cctype>
#include <ostream>
#include <string>
#include <vector>

namespace gfe::experiment {

/**
 * Compute a set of basic statistics (mean, median, and so on) from a sequence of results.
 * The class is not thread safe
 */
class ExecStatistics {
    friend std::ostream& operator<<(std::ostream& out, const ExecStatistics& stats);

    const uint64_t m_num_trials;
    uint64_t m_num_timeouts {0};
    uint64_t m_mean {0};
    uint64_t m_stddev {0};
    uint64_t m_min {0};
    uint64_t m_max {0};
    uint64_t m_median {0};
    uint64_t m_percentile90 {0};
    uint64_t m_percentile95 {0};
    uint64_t m_percentile97 {0};
    uint64_t m_percentile99 {0};

private:
    static uint64_t get_percentile(const std::vector<uint64_t>& values_sorted, uint64_t position);

public:
    /**
     * Compute the statistics for the given sequence of results. Each entry represents either
     * the completion time (in microsecs), or the special value -1 to indicate a timeout.
     */
    ExecStatistics(const std::vector<int64_t>& trials);

    /**
     * Save the computed statistics in the database
     */
    void save(const std::string& name);
};

// Print the statistics into the given output stream, for reporting or debugging purposes
std::ostream& operator<<(std::ostream& out, const ExecStatistics& stats);

} // namespace



