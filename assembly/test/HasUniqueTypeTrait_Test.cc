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