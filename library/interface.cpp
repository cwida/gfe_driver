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

#include "interface.hpp"

#include <memory>

#if defined(HAVE_STINGER)
#include "stinger/stinger.hpp"
#endif

using namespace std;

namespace library {

/*****************************************************************************
 *                                                                           *
 *  Factory                                                                  *
 *                                                                           *
 *****************************************************************************/
ImplementationManifest::ImplementationManifest(const string& name, const string& description, unique_ptr<Interface> (*factory)(void)) :
    m_name(name), m_description(description), m_factory(factory){ }

#if defined(HAVE_STINGER)
std::unique_ptr<Interface> generate_stinger(){
    return unique_ptr<Interface>{ new Stinger() };
}
#endif

vector<ImplementationManifest> implementations() {
    vector<ImplementationManifest> result;

#if defined(HAVE_STINGER)
    result.emplace_back("stinger", "Stinger library", &generate_stinger);
#endif

    return result;
}

/*****************************************************************************
 *                                                                           *
 *  Base interface                                                           *
 *                                                                           *
 *****************************************************************************/
Interface::Interface(){}
Interface::~Interface(){}
void Interface::on_main_init(int num_threads){ };
void Interface::on_thread_init(int thread_id){ };
void Interface::on_thread_destroy(int thread_id){ } ;
void Interface::on_main_destroy(){ };

} // namespace library
