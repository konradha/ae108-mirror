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

#include "ae108/assembly/FeaturePlugins.h"

namespace ae108 {
namespace assembly {

/**
 * @brief ConcatenatePlugins is a metaprogramming function that concatenates
 * FeaturePlugins. It contains a "type" typedef of the concatenated
 * FeaturePlugins.
 */
template <class... FeaturePlugins> struct ConcatenatePlugins {};

} // namespace assembly
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

namespace ae108 {
namespace assembly {
namespace detail {

/**
 * @brief A helper class for plugin list concatenation.
 */
template <class, class> struct ConcatenatePluginsImpl {};

/**
 * @brief Concatenates exactly two lists.
 */
template <template <class> class... Lhs, template <class> class... Rhs>
struct ConcatenatePluginsImpl<FeaturePlugins<Lhs...>, FeaturePlugins<Rhs...>> {
  using type = FeaturePlugins<Lhs..., Rhs...>;
};
} // namespace detail

/**
 * @brief Concatenate an empty list.
 */
template <> struct ConcatenatePlugins<> { using type = FeaturePlugins<>; };

/**
 * @brief Concatenate a single list.
 */
template <template <class> class... Plugins>
struct ConcatenatePlugins<FeaturePlugins<Plugins...>> {
  using type = FeaturePlugins<Plugins...>;
};

/**
 * @brief Concatenate two or more lists.
 */
template <class List, class... Lists>
struct ConcatenatePlugins<List, Lists...> {
  using type = typename detail::ConcatenatePluginsImpl<
      List, typename ConcatenatePlugins<Lists...>::type>::type;
};
} // namespace assembly
} // namespace ae108
