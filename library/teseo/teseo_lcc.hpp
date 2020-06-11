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

#include "teseo_driver.hpp"

namespace gfe::library {

/**
 * Driver for Teseo, with a specialised implementation of the LCC kernel
 */
class TeseoLCC : public TeseoDriver {
    TeseoLCC(const TeseoLCC& ) = delete;
    TeseoLCC& operator=(const TeseoLCC& ) = delete;
    const bool m_low_level_api; // use the low level API to implement LCC?

public:
    /**
     * Constructor, same arguments of its base class
     */
    TeseoLCC(bool is_directed, bool low_level_api);

    /**
     * Specialised implementation of the kernel LCC
     */
    virtual void lcc(const char* dump2file = nullptr);
};

} // namespace
