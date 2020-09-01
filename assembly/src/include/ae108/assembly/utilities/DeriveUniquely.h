// © 2020 ETH Zurich, Mechanics and Materials Lab
// © 2020 California Institute of Technology
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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