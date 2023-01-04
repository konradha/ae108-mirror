// © 2020 ETH Zurich, Mechanics and Materials Lab
// © 2020 California Institute of Technology
//
// This file is part of ae108.
//
// ae108 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or any
// later version.
//
// ae108 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ae108. If not, see <https://www.gnu.org/licenses/>.

#pragma once

#include <type_traits>

namespace ae108 {
namespace assembly {
namespace utilities {

/**
 * @brief Inherits std::true_type if Args contains exactly one type,
 * disregarding duplicates.
 */
template <class... Args> struct IsSingleUniqueType : std::true_type {};

template <> struct IsSingleUniqueType<> : std::false_type {};

template <class A> struct IsSingleUniqueType<A> : std::true_type {
  using type = A;
};

template <class A, class B, class... Args>
struct IsSingleUniqueType<A, B, Args...>
    : std::conditional<std::is_same<A, B>::value,
                       IsSingleUniqueType<A, Args...>, std::false_type>::type {
};
} // namespace utilities
} // namespace assembly
} // namespace ae108