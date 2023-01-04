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
