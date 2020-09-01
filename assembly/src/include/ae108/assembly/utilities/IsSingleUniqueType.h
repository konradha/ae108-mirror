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