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

namespace ae108 {
namespace assembly {
namespace utilities {

/**
 * @brief Derives each of its template arguments once.
 */
template <class... Args> struct DeriveUniquely;
} // namespace utilities
} // namespace assembly
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

#include <type_traits>

namespace ae108 {
namespace assembly {
namespace utilities {
namespace detail {
template <class... Tags> struct Empty {};

template <class Base, class Derived, class... Tags>
using EmptyIfAlreadyBase =
    std::conditional<std::is_base_of<Base, Derived>::value, Empty<Tags...>,
                     Base>;
} // namespace detail

template <> struct DeriveUniquely<> {};

template <class Arg, class... Args>
struct DeriveUniquely<Arg, Args...>
    : DeriveUniquely<Args...>,
      detail::EmptyIfAlreadyBase<Arg, DeriveUniquely<Args...>, Args...>::type {
};
} // namespace utilities
} // namespace assembly
} // namespace ae108