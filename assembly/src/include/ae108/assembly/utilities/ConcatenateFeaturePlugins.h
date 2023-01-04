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
