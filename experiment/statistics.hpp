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
    uint64_t m_mean {0};
    uint64_t m_stddev {0};
    uint64_t m_min {0};
    uint64_t m_max {0};
    uint64_t m_percentile[99]; // percentile 1, ..., 99

public:
    ExecStatistics(std::vector<uint64_t>& trials);

    /**
     * Save in the database
     */
    void save(const std::string& name);
};

std::ostream& operator<<(std::ostream& out, const ExecStatistics& stats);

} // namespace



