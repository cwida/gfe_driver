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

/**
 * The library requires to be set up with a few macros to enable/disable its features. Which is
 * what we do here.
 * This header acts a wrapper before including the actual llama.h header, ensuring all translation
 * units are compiled in a consistent way.
 */

// Which implementation to use?
#define LL_MEMORY_ONLY /* llama, ask for the implementation in-memory */
//#define LL_PERSISTENT /* llama, ask for the implementation with memory mapped files */

// Support for deletions?
#define LL_DELETIONS

// Disable a warning in the llama.h implementation, ISO C++17 does not allow register bla bla
#define register /* nop */

// We're finally ready to include the actual library
#include "llama.h"

// the string ID for the vertex property associated to the external IDs (e.g. user id)
extern char const * const g_llama_property_names;

// the string ID for the edge property associated to the weights
extern char const * const g_llama_property_weights;
