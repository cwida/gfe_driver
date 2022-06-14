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

#include "quantity.hpp"

#include <cassert>
#include <cctype> // tolower
#include <cmath>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>

#undef CURRENT_ERROR_TYPE
#define CURRENT_ERROR_TYPE QuantityError

using namespace common;
using namespace std;

/**********************************************************************************************************************
 *                                                                                                                    *
 *  ComputerQuantity                                                                                                  *
 *                                                                                                                    *
 **********************************************************************************************************************/

ComputerQuantity::ComputerQuantity() : ComputerQuantity(0) {}

ComputerQuantity::ComputerQuantity(uint64_t magnitude, bool byte_suffix) : m_magnitude(magnitude), m_is_byte_quantity(byte_suffix) { }

ComputerQuantity::ComputerQuantity(const string& quantity, bool byte_suffix) :
        m_magnitude(parse(quantity, byte_suffix ? ByteSuffix::Optional : ByteSuffix::Missing)),
        m_is_byte_quantity(byte_suffix){ }

uint64_t ComputerQuantity::magnitude() const noexcept {
    return m_magnitude;
}

ComputerQuantity::operator int64_t() const noexcept {
    return static_cast<int64_t>(m_magnitude);
}

uint64_t ComputerQuantity::parse(const string& quantity, ByteSuffix byte_suffix) {
    string pattern {"^\\s*(\\d+)?(\\.\\d+)?\\s*([KMGT]?)("}; // kilo, mega, giga tera
    if(byte_suffix == ByteSuffix::Optional || byte_suffix == ByteSuffix::Mandatory) pattern += "B"; // bytes
    if(byte_suffix == ByteSuffix::Optional) pattern += "?";
    pattern += ")\\s*$";
    regex regex{pattern, regex::icase};
    smatch match;
    if (!regex_match(begin(quantity), end(quantity), match, regex)){ ERROR("Invalid quantity: " << quantity); }

    // check the `byte' suffix is present
    if(byte_suffix == ByteSuffix::Mandatory && match[4].length() == 0){
        ERROR("Invalid quantity: " << quantity << ". The byte (b) suffix is required.");
    }

    // multiply the quantity by the unit
    uint64_t mult = 1;
    if(match[3].length() > 0){
        assert(match[3].length() == 1 && "Expected one of the following single units K, M, G or T");
        char unit = tolower(match[3].str().c_str()[0]);
        switch(unit){
        case 'k': mult = (1ull <<10); break;
        case 'm': mult = (1ull <<20); break;
        case 'g': mult = (1ull <<30); break;
        case 't': mult = (1ull <<40); break;
        default: assert(0 && "Invalid unit");
        }
    }

    // parse the magnitude
    uint64_t result{0};
    string str_value = match[1];
    str_value += match[2];
    if(str_value.empty()){ ERROR("Magnitude missing: " << quantity); }
    bool is_decimal = match[2].length() > 0;
    if(is_decimal){
        size_t leftover;
        double value = stod(str_value, &leftover);
        if(leftover != str_value.length()) ERROR("Parse error: " << quantity);
        value *= mult;
        double intpart;
        double fractpart = modf(value, &intpart);
        if(fractpart != 0.) ERROR("Cannot cast to an integer value: " << quantity << ". Absolute value: " << value);
        result = value;
    } else {
        size_t leftover;
        result = stoull(str_value, &leftover);
        if(leftover != str_value.length()) ERROR("Parse error: " << quantity);
        result *= mult;
    }

    return result;
}

static string convert_to_string(uint64_t magnitude, ComputerQuantity::Unit unit, bool byte_suffix){
    using Unit = ComputerQuantity::Unit;

    // unit
    uint64_t mult = 1;
    string suffix;
    switch(unit){
    case Unit::BASIC: break; // nop
    case Unit::KILO: mult = (1ull << 10); suffix = "K"; break;
    case Unit::MEGA: mult = (1ull << 20); suffix = "M"; break;
    case Unit::GIGA: mult = (1ull << 30); suffix = "G"; break;
    case Unit::TERA: mult = (1ull << 40); suffix = "T"; break;
    default: ERROR("Invalid unit: " << (int) unit);
    }
    if(byte_suffix){
        if(mult == 1){
            suffix += "bytes";
        } else {
            suffix += "B";
        }
    }

    stringstream ss;
    if(magnitude % mult == 0){
        ss << magnitude / mult;
    } else {
        char buffer[128];
        snprintf(buffer, 128, "%.2f", ((double) magnitude) / mult);
        ss << buffer;
    }

    if(!suffix.empty()) {
        if (byte_suffix) ss << " ";
        ss << suffix;
    }

    return ss.str();
}

string ComputerQuantity::to_string(ComputerQuantity::Unit unit) const {
    switch(unit){
    case Unit::AUTO:
        if(m_magnitude >= (1ull << 40)){
            return to_string(Unit::TERA);
        } else if(m_magnitude >= (1ull << 30)){
            return to_string(Unit::GIGA);
        } else if(m_magnitude >= (1ull << 20)){
            return to_string(Unit::MEGA);
        } else if(m_magnitude >= (1ull << 10)){
            return to_string(Unit::KILO);
        } else {
            return to_string(Unit::BASIC);
        }
    default:
        return convert_to_string(m_magnitude, unit, m_is_byte_quantity);
    }
}

ostream& (::common::operator<<)(ostream& out, const ComputerQuantity& q) {
    out << q.to_string();
    return out;
}

ostream& (::common::operator<<)(ostream& out, const ComputerQuantity* q){
    out << "ptr: 0x" << std::hex << (void*) q << std::dec;
    if(q != nullptr){
        out << " (" << q->to_string() << ")";
    }
    return out;
}

std::istream& (::common::operator>>)(std::istream& in, common::ComputerQuantity& q){
    string value;
    in >> value;
    auto magnitude = ComputerQuantity::parse(value, q.is_byte_quantity() ? ComputerQuantity::ByteSuffix::Mandatory : ComputerQuantity::ByteSuffix::Missing);
    q.m_magnitude = magnitude;
    return in;
}

ComputerQuantity (::common::operator+)(ComputerQuantity q1, int64_t q2) {
    int64_t v1 = q1.magnitude();
    int64_t v2 = q2;
    int64_t res = v1 + v2;
    if(res < 0) ERROR("Negative quantity: " << v1 << " + " << v2);
    return ComputerQuantity(res, q1.is_byte_quantity());
}

ComputerQuantity (::common::operator-)(ComputerQuantity q1, int64_t q2){
    return q1 + (-q2);
}

ComputerQuantity (::common::operator*)(ComputerQuantity q1, int64_t q2){
    int64_t v1 = q1.magnitude();
    int64_t v2 = q2;
    int64_t res = v1 * v2;
    if(res < 0) ERROR("Negative quantity: " << q1 << " * " << q2);
    return ComputerQuantity(res, q1.is_byte_quantity());
}

ComputerQuantity (::common::operator/)(ComputerQuantity q1, int64_t q2){
    int64_t n = q1.magnitude();
    int64_t d = q2;
    int64_t res = n / d;
    if(res < 0) ERROR("Negative quantity: " << q1 << " * " << q2);
    return ComputerQuantity(res, q1.is_byte_quantity());
}

ComputerQuantity& ComputerQuantity::operator+=(int64_t value) {
    int64_t v1 = magnitude();
    int64_t v2 = value;
    int64_t res = v1 + v2;
    if(res < 0) ERROR("Negative quantity: " << v1 << " + " << v2);
    m_magnitude = static_cast<uint64_t>(res);
    return *this;
}

ComputerQuantity& ComputerQuantity::operator-=(int64_t value) {
    return this->operator+=(-value);
}

ComputerQuantity& ComputerQuantity::operator*=(int64_t v2) {
    int64_t v1 = magnitude();
    int64_t res = v1 * v2;
    if(res < 0) ERROR("Negative quantity: " << *this << " * " << v2);
    m_magnitude = static_cast<uint64_t>(res);
    return *this;
}

ComputerQuantity& ComputerQuantity::operator/=(int64_t d) {
    int64_t n = magnitude();
    int64_t res = n / d;
    if(res < 0) ERROR("Negative quantity: " << *this << " * " << d);
    m_magnitude = static_cast<uint64_t>(res);
    return *this;
}

/**********************************************************************************************************************
 *                                                                                                                    *
 *  DurationQuantity                                                                                                  *
 *                                                                                                                    *
 **********************************************************************************************************************/
DurationQuantity::DurationQuantity() : m_duration(0){
    /* nop */
}

DurationQuantity::DurationQuantity(uint64_t seconds) : DurationQuantity(chrono::seconds(seconds)) {
    /* nop */
}

DurationQuantity::operator uint64_t() const noexcept {
    return seconds();
}

uint64_t DurationQuantity::seconds() const noexcept {
    return chrono::duration_cast<chrono::seconds>(m_duration).count();
}

DurationQuantity DurationQuantity::parse(const std::string& str_duration){
    static string pattern { R"raw(^(\d+)?\s*(\w+)?\s*$)raw" }; // kilo, mega, giga tera
    enum class Unit { HOURS, MINUTES, SECONDS, MILLISECONDS, MICROSECONDS, NANOSECONDS };
    struct Association { Unit m_unit; string m_string; };
    static Association associations[] = {
            {Unit::HOURS, "h"},
            {Unit::HOURS, "hour"},
            {Unit::HOURS, "hours"},
            {Unit::MINUTES, "m"},
            {Unit::MINUTES, "min"},
            {Unit::MINUTES, "mins"},
            {Unit::MINUTES, "minute"},
            {Unit::MINUTES, "minutes"},
            {Unit::SECONDS, "s"},
            {Unit::SECONDS, "sec"},
            {Unit::SECONDS, "secs"},
            {Unit::SECONDS, "second"},
            {Unit::SECONDS, "seconds"},
            {Unit::MILLISECONDS, "ms"},
            {Unit::MILLISECONDS, "msec"},
            {Unit::MILLISECONDS, "msecs"},
            {Unit::MILLISECONDS, "millisec"},
            {Unit::MILLISECONDS, "millisecs"},
            {Unit::MILLISECONDS, "millisecond"},
            {Unit::MILLISECONDS, "milliseconds"},
            {Unit::MICROSECONDS, "us"},
            {Unit::MICROSECONDS, "usec"},
            {Unit::MICROSECONDS, "usecs"},
            {Unit::MICROSECONDS, "microsec"},
            {Unit::MICROSECONDS, "microsecs"},
            {Unit::MICROSECONDS, "microsecond"},
            {Unit::MICROSECONDS, "microseconds"},
            {Unit::NANOSECONDS, "ns"},
            {Unit::NANOSECONDS, "nsec"},
            {Unit::NANOSECONDS, "nsecs"},
            {Unit::NANOSECONDS, "nanosec"},
            {Unit::NANOSECONDS, "nanosecs"},
            {Unit::NANOSECONDS, "nanosecond"},
            {Unit::NANOSECONDS, "nanoseconds"},
    };

    regex regex{pattern};
    smatch match;
    if (!regex_match(begin(str_duration), end(str_duration), match, regex)){ ERROR("Invalid duration: " << str_duration); }

    // if a unit has not been explicitly specified, assume seconds
    Unit unit = Unit::SECONDS;
    int64_t duration = stoi( match[1] );
    if(duration < 0) ERROR("Invalid duration: the given value is negative: " << duration);
    if(match[2].length() > 0){
        string str_unit = match[2];
        transform(begin(str_unit), end(str_unit), begin(str_unit), ::tolower); // lower case
        uint64_t num_associations = sizeof(associations) / sizeof(associations[0]);
        bool stop = false;
        uint64_t i = 0;
        while(i < num_associations && !stop){
            if(str_unit == associations[i].m_string){
                unit = associations[i].m_unit;
                stop = true;
            } else {
                i++;
            }
        }
        if(!stop){
            ERROR("Invalid duration unit: `" << match[2] << "'");
        }
    }

    switch(unit){
    case Unit::NANOSECONDS:
        return DurationQuantity(chrono::nanoseconds(duration));
    case Unit::MICROSECONDS:
        return DurationQuantity(chrono::microseconds(duration));
    case Unit::MILLISECONDS:
        return DurationQuantity(chrono::milliseconds(duration));
    case Unit::SECONDS:
        return DurationQuantity(chrono::seconds(duration));
    case Unit::MINUTES:
        return DurationQuantity(chrono::minutes(duration));
    case Unit::HOURS:
        return DurationQuantity(chrono::hours(duration));
    default:
        ERROR("Unit not handled: " << (int) unit);
    }
}

std::string DurationQuantity::to_string() const {
    uint64_t nanoseconds = as<chrono::nanoseconds>().count();
    uint64_t microseconds = nanoseconds / 1000;
    if(microseconds > 0) nanoseconds %= 1000;
    uint64_t milliseconds = microseconds/ 1000;
    if(milliseconds > 0) microseconds %= 1000;
    uint64_t seconds = milliseconds / 1000;
    if(seconds > 0) milliseconds %= 1000;
    uint64_t minutes = seconds / 60;
    if(minutes > 0) seconds %= 60;
    uint64_t hours = minutes / 60;
    if(hours > 0) minutes %= 60;

    stringstream ss;

    if(hours > 0) {
        ss << hours << " hour";
        if(hours > 1){ ss << "s"; }

        if (minutes > 0) {
            if (seconds > 0) {
                ss << ", ";
            } else {
                ss << " and ";
            }
            ss << minutes << " minute";
            if(minutes > 1) { ss << "s"; }
        }
        if (seconds > 0) {
            ss << " and " << seconds << " second";
            if(seconds > 1){ ss << "s"; }
        }
    } else if (minutes > 0){
        ss << minutes << " minute";
        if(minutes > 1) { ss << "s"; }
        if(seconds > 0){
            ss << " and " << seconds << " second";
            if(seconds > 1){ ss << "s"; }
        }
    } else if (seconds > 0){
        ss << seconds;
        if(milliseconds > 0){
            ss << "." << milliseconds;
        }
        ss << " second";
        if (seconds > 1) { ss << "s"; }
    } else if (milliseconds > 0){
        ss << milliseconds;
        if(milliseconds < 10 && microseconds > 0){
            ss << "." << microseconds;
        }
        ss << " millisecs";
    } else if (microseconds > 0){
        ss << microseconds;
        if(microseconds < 10 && nanoseconds > 0){
            ss << "." << nanoseconds;
        }
        ss << " microsecs";
    } else { // nanoseconds
        ss << nanoseconds << " nanoseconds";
    }

    return ss.str();
}

ostream& (::common::operator<<)(ostream& out, const DurationQuantity& q) {
    out << q.to_string();
    return out;
}

ostream& (::common::operator<<)(ostream& out, const DurationQuantity* q){
    out << "ptr: 0x" << std::hex << (void*) q << std::dec;
    if(q != nullptr){
        out << " (" << q->to_string() << ")";
    }
    return out;
}

std::istream& (::common::operator>>)(std::istream& in, common::DurationQuantity& q){
    string value;
    getline(in, value);
    auto duration_parsed = DurationQuantity::parse(value);
    q.m_duration = duration_parsed.m_duration;
    return in;
}
