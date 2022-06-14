/**
 * This file is part of libcommon.
 *
 * libcommon is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, orF
 * (at your option) any later version.
 *
 * libcommon is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with libcommon.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "system.hpp"

#include <iostream>
#include <sstream>
#include <string>

using namespace std;

namespace common {

CompilerFamily CompilerInfo::family() const {
    return m_family;
}

string CompilerInfo::name() const {
    return m_name;
}

int CompilerInfo::major() const {
    return m_major;
}

int CompilerInfo::minor() const {
    return m_minor;
}

int CompilerInfo::patch() const {
    return m_patch;
}

string CompilerInfo::version() const {
    stringstream ss;
    ss << m_major << "." << m_minor;
    if(m_patch != 0){ ss << "." << m_patch; }
    return ss.str();
}

string CompilerInfo::to_string() const {
    stringstream ss;

    if(m_family == CompilerFamily::UNKNOWN){
        ss << "[compiler unknown]";
    } else {
        ss << name() << " " << version();
    }

    return ss.str();
}

std::ostream& operator<<(std::ostream& out, const CompilerInfo& ci){
    out << ci.to_string();
    return out;
}

} // namespace



