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

namespace experiment {

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
    ExecStatistics(const std::vector<int64_t>& trials);

    /**
     * Save in the database
     */
    void save(const std::string& name);
};

std::ostream& operator<<(std::ostream& out, const ExecStatistics& stats);

} // namespace



