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
