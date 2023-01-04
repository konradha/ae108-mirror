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

#include "ae108/assembly/AssemblerTypeTraits.h"
#include "ae108/assembly/utilities/ConcatenateFeaturePlugins.h"
#include "ae108/assembly/utilities/DerivePluginsUniquely.h"
#include "ae108/assembly/utilities/HasUniqueTypeTrait.h"
#include <cstddef>
#include <tuple>
#include <type_traits>
#include <utility>

namespace ae108 {
namespace assembly {

/**
 * @brief Groups a list of single element assemblers and inherits all their
 * plugins. When calling a method provided by a plugin, this method is called
 * for all member assemblers that have that plugin.
 */
template <class... SingleElementAssemblers>
class AssemblerGroup : public DerivePluginsUniquely<
                           AssemblerGroup<SingleElementAssemblers...>,
                           typename ConcatenatePlugins<typename PluginTypeTrait<
                               SingleElementAssemblers>::type...>::type> {
public:
  template <class... Args> explicit AssemblerGroup(Args &&...args);

  /**
   * @brief Get the Nth member assembler.
   */
  template <std::size_t N> const auto &get() const & {
    return std::get<N>(_assemblers);
  }

  /**
   * @brief Get the Nth member assembler.
   */
  template <std::size_t N> auto &get() & { return std::get<N>(_assemblers); }

  /**
   * @brief Return the tuple of assemblers.
   */
  const auto &assemblers() const &noexcept { return _assemblers; }

  /**
   * @brief Return the tuple of assemblers.
   */
  auto &assemblers() &noexcept { return _assemblers; }

private:
  using tuple_type = std::tuple<SingleElementAssemblers...>;
  tuple_type _assemblers;
};

/**
 * @brief Specialization: Deduces the policy type corresponding to an Assembler
 * group.
 */
template <class... SingleElementAssemblers>
struct PolicyTypeTrait<AssemblerGroup<SingleElementAssemblers...>>
    : utilities::HasUniqueTypeTrait<PolicyTypeTrait,
                                    SingleElementAssemblers...> {};

/**
 * @brief Specialization: Deduces the mesh type corresponding to an Assembler
 * group.
 */
template <class... SingleElementAssemblers>
struct MeshTypeTrait<AssemblerGroup<SingleElementAssemblers...>>
    : utilities::HasUniqueTypeTrait<MeshTypeTrait, SingleElementAssemblers...> {
};

/**
 * @brief Specialization: Deduces the vector type corresponding to an Assembler
 * group.
 */
template <class... SingleElementAssemblers>
struct VectorTypeTrait<AssemblerGroup<SingleElementAssemblers...>>
    : utilities::HasUniqueTypeTrait<VectorTypeTrait,
                                    SingleElementAssemblers...> {};

/**
 * @brief Specialization: Deduces the matrix type corresponding to an Assembler
 * group.
 */
template <class... SingleElementAssemblers>
struct MatrixTypeTrait<AssemblerGroup<SingleElementAssemblers...>>
    : utilities::HasUniqueTypeTrait<MatrixTypeTrait,
                                    SingleElementAssemblers...> {};

/**
 * @brief Specialization: Deduces the size type corresponding to an Assembler
 * group.
 */
template <class... SingleElementAssemblers>
struct SizeTypeTrait<AssemblerGroup<SingleElementAssemblers...>>
    : utilities::HasUniqueTypeTrait<SizeTypeTrait, SingleElementAssemblers...> {
};

/**
 * @brief Specialization: Deduces the value type corresponding to an Assembler
 * group.
 */
template <class... SingleElementAssemblers>
struct ValueTypeTrait<AssemblerGroup<SingleElementAssemblers...>>
    : utilities::HasUniqueTypeTrait<ValueTypeTrait,
                                    SingleElementAssemblers...> {};

/**
 * @brief Specialization: Deduces the value type corresponding to an Assembler
 * group.
 */
template <class... SingleElementAssemblers>
struct RealTypeTrait<AssemblerGroup<SingleElementAssemblers...>>
    : utilities::HasUniqueTypeTrait<ValueTypeTrait,
                                    SingleElementAssemblers...> {};

/**
 * @brief Specialization: Deduces the element type corresponding to an Assembler
 * group.
 */
template <class... SingleElementAssemblers>
struct ElementTypeTrait<AssemblerGroup<SingleElementAssemblers...>>
    : utilities::HasUniqueTypeTrait<ElementTypeTrait,
                                    SingleElementAssemblers...> {};

/**
 * @brief Specialization: AssemblerGroup is an assembler group.
 */
template <class... SingleElementAssemblers>
struct IsGroupTypeTrait<AssemblerGroup<SingleElementAssemblers...>>
    : std::true_type {};

} // namespace assembly
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

namespace ae108 {
namespace assembly {

template <class... SingleElementAssemblers>
template <class... Args>
AssemblerGroup<SingleElementAssemblers...>::AssemblerGroup(Args &&...args)
    : _assemblers(std::forward<Args>(args)...) {}

} // namespace assembly
} // namespace ae108
