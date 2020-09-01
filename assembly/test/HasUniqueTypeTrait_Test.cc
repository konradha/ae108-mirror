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

#include "ae108/assembly/utilities/HasUniqueTypeTrait.h"
#include <type_traits>

namespace ae108 {
namespace assembly {
namespace utilities {
namespace {

struct A {};
struct B {};
struct C {};
struct D {};

template <class T> struct trait {};

template <> struct trait<A> { using type = C; };
template <> struct trait<B> { using type = D; };

static_assert(HasUniqueTypeTrait<trait, A>::value, "C is single unique type.");

static_assert(!HasUniqueTypeTrait<trait>::value,
              "Empty is not a single unique type.");

static_assert(HasUniqueTypeTrait<trait, A, A>::value,
              "C, C is a single unique type.");

static_assert(
    std::is_same<HasUniqueTypeTrait<trait, A, A>::type, trait<A>::type>::value,
    "HasUniqueTypeTrait returns trait<A>.");

static_assert(!HasUniqueTypeTrait<trait, A, A, B>::value,
              "C, C, D is not a single unique type.");

static_assert(!HasUniqueTypeTrait<trait, A, B>::value,
              "C, D is not a single unique type.");

static_assert(!HasUniqueTypeTrait<trait, A, B, A>::value,
              "C, D, C is not a single unique type.");

} // namespace
} // namespace utilities
} // namespace assembly
} // namespace ae108