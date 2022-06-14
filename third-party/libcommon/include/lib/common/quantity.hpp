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

#ifndef COMMON_QUANTITY_HPP
#define COMMON_QUANTITY_HPP

#include <cinttypes>
#include <chrono>
#include <cstddef>
#include <iostream>
#include <string>

#include "error.hpp"

namespace common {

/**
 * The kind of exception thrown by the members of the class Quantity
 */
DEFINE_EXCEPTION(QuantityError);

/**
 * It represents a computer quantity, e.g. 4 K or 8 GB, where 1k = 1024
 */
class ComputerQuantity {
    friend ::std::istream& operator>>(::std::istream& in, ::common::ComputerQuantity& q);
private:
    uint64_t m_magnitude; // the actual magnitude, either in bytes or as absolute number
    bool m_is_byte_quantity; // append the bytes suffix

public:
    /**
     * Empty quantity (zero)
     */
    ComputerQuantity();

    ComputerQuantity(uint64_t magnitude, bool byte_suffix = false);

    ComputerQuantity(const std::string& quantity, bool byte_suffix = false);

    /**
     * The actual magnitude of this quantity, in its most basic unit (e.g. bytes/absolute number).
     */
    uint64_t magnitude() const noexcept;

    /**
     * Return the absolute magnitude, in its basic unit, associated to this quantity.
     * Cast to int, it is usually simpler to handle.
     */
    operator int64_t() const noexcept;

    /**
     * Whether this is a byte quantity, that is has the suffix "B" for bytes, or an absolute number, e.g. 2 K
     */
    bool is_byte_quantity() const {
        return m_is_byte_quantity;
    }

    /**
     * Set the flag for the byte suffix (B)
     */
    void set_byte_quantity(bool value) {
        m_is_byte_quantity = value;
    }

    /**
     * Whether the string to parse can or should have the suffix 'b' for bytes
     */
    enum class ByteSuffix { Missing, Optional, Mandatory };

    /**
     * Parse the given string and return the absolute magnitude in its basic unit
     */
    static uint64_t parse(const std::string& quantity, ByteSuffix byte_suffix);


    enum class Unit { AUTO, BASIC, KILO, MEGA, GIGA, TERA };

    /**
     * Return a string representation of the quantity
     * @param unit the unit to use
     * @return a string representation of the quantity
     */
    std::string to_string(Unit unit = Unit::AUTO) const;

    /**
     * Mathematical operations
     */
    ComputerQuantity& operator+=(int64_t value);
    ComputerQuantity& operator-=(int64_t value);
    ComputerQuantity& operator*=(int64_t value);
    ComputerQuantity& operator/=(int64_t value);
};

ComputerQuantity operator+(ComputerQuantity q1, int64_t q2);
ComputerQuantity operator-(ComputerQuantity q1, int64_t q2);
ComputerQuantity operator*(ComputerQuantity q1, int64_t q2);
ComputerQuantity operator/(ComputerQuantity q1, int64_t q2);

::std::ostream& operator<<(::std::ostream& out, const ComputerQuantity& q);
::std::ostream& operator<<(::std::ostream& out, const ComputerQuantity* q);
::std::istream& operator>>(::std::istream& in, ComputerQuantity& q);

/**
* It represents a time quantity, e.g. 4h, 8 secs or 30 microsecs
*/
class DurationQuantity {
    friend ::std::istream& operator>>(::std::istream& in, ::common::DurationQuantity& q);
    std::chrono::nanoseconds m_duration; // internal representation

public:
    /**
     * Empty quantity (zero)
     */
    DurationQuantity();

    /**
     * Init the quantity in terms of seconds
     */
     explicit DurationQuantity(uint64_t seconds);

    /**
     * Init the quantity with the given duration
     */
     template<typename D>
     explicit DurationQuantity(D duration);

    /**
     * Return the duration in seconds
     */
    operator uint64_t() const noexcept;

    /**
     * Return the duration in seconds
     */
    uint64_t seconds() const noexcept;

    /**
     * Parse the given string and return the absolute magnitude in its basic unit
     */
    static DurationQuantity parse(const std::string& str_duration);

    /**
     * Return a string representation of the duration
     */
    std::string to_string() const;

    /**
     * Return the time with the given representation D
     */
     template<typename D>
     D as() const;
};

::std::ostream& operator<<(::std::ostream& out, const DurationQuantity& q);
::std::ostream& operator<<(::std::ostream& out, const DurationQuantity* q);
::std::istream& operator>>(::std::istream& in, DurationQuantity& q);

// implementation details
template<typename D>
DurationQuantity::DurationQuantity(D duration) : m_duration(std::chrono::duration_cast<std::chrono::nanoseconds>(duration)) {
    /* nop */
}

template<typename D>
D DurationQuantity::as() const {
    return std::chrono::duration_cast<D>(m_duration);
}

} // namespace common

#endif //COMMON_QUANTITY_HPP
