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

#include <cinttypes>
#include <memory>

namespace gfe::graph { class WeightedEdgeStream; } // forward declaration
namespace gfe::library { class Interface; } // forward declaration

namespace gfe::experiment {

/**
 * Check that all edges in the stream are contained in the interface. Report the number of missing vertices (0 => validation successful).
 */
uint64_t validate_updates(std::shared_ptr<gfe::library::Interface> interface, std::shared_ptr<gfe::graph::WeightedEdgeStream> stream);

} // namespace
