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

#include "ae108/assembly/FeaturePlugins.h"
#include "ae108/assembly/utilities/ConcatenateFeaturePlugins.h"
#include <type_traits>

namespace ae108 {
namespace assembly {
namespace {

template <class> struct A {};
template <class> struct B {};
template <class> struct C {};

static_assert(std::is_same<ConcatenatePlugins<>::type, FeaturePlugins<>>::value,
              "Concatenating no list yields empty list.");

static_assert(std::is_same<ConcatenatePlugins<FeaturePlugins<A, B>>::type,
                           FeaturePlugins<A, B>>::value,
              "Concatenating single list yields same list.");

static_assert(std::is_same<ConcatenatePlugins<FeaturePlugins<A, B>,
                                              FeaturePlugins<B, C>>::type,
                           FeaturePlugins<A, B, B, C>>::value,
              "Concatenating two lists yields concatenated list.");

static_assert(
    std::is_same<ConcatenatePlugins<FeaturePlugins<A>, FeaturePlugins<B>,
                                    FeaturePlugins<C>>::type,
                 FeaturePlugins<A, B, C>>::value,
    "Concatenating three lists yields concatenated list.");

} // namespace
} // namespace assembly
} // namespace ae108
