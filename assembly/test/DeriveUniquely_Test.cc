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

#include "ae108/assembly/utilities/DeriveUniquely.h"
#include <type_traits>

namespace ae108 {
namespace assembly {
namespace utilities {
namespace {

struct A {};
struct B {};

struct C : DeriveUniquely<A, B, A> {};

static_assert(std::is_base_of<A, C>::value, "A is a base of C.");

static_assert(std::is_base_of<B, C>::value, "B is a base of C.");

struct EmptyListCompiles : DeriveUniquely<> {};

} // namespace
} // namespace utilities
} // namespace assembly
} // namespace ae108