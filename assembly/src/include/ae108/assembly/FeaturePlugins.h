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

#include <utility>

namespace ae108 {
namespace assembly {

namespace detail {
template <class Assembler, template <class> class... Plugins>
struct FeaturePluginImpl : Plugins<Assembler>... {
  /**
   * @brief Pass no variable to every plugin.
   */
  explicit FeaturePluginImpl() = default;

  /**
   * @brief Copy-initialize every plugin.
   */
  explicit FeaturePluginImpl(Plugins<Assembler>... plugins)
      : Plugins<Assembler>(std::move<Plugins<Assembler>>(plugins))... {}

protected:
  ~FeaturePluginImpl() = default;
};

/**
 * @brief Specialization for no plugins.
 */
template <class Assembler> struct Empty {
protected:
  ~Empty() = default;
};
} // namespace detail

/**
 * @brief A wrapper around a list of feature plugins.
 */
template <template <class> class... Plugins> struct FeaturePlugins {
  template <class Assembler>
  using type = detail::FeaturePluginImpl<Assembler, Plugins...>;
};

/**
 * @brief Specialization for no plugins.
 */
template <> struct FeaturePlugins<> {
  template <class Assembler> using type = detail::Empty<Assembler>;
};
} // namespace assembly
} // namespace ae108
