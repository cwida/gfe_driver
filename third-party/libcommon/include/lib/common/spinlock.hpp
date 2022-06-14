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

#ifndef COMMON_SPINLOCK_HPP
#define COMMON_SPINLOCK_HPP

#include <atomic>

namespace common {

class SpinLock {
    std::atomic_flag m_lock = ATOMIC_FLAG_INIT; // atomic integral
public:
    /**
     * Acquire the lock
     */
    void lock() {
        while (m_lock.test_and_set(std::memory_order_acquire)) /* nop */;
    }

    /**
     * Release the lock
     */
    void unlock() {
        m_lock.clear(std::memory_order_release);
    }
};

} // namespace common

#endif /* COMMON_SPINLOCK_HPP */