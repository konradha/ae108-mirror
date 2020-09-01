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
  template <std::size_t N>
  using member_assembler_type =
      typename MemberTypeTrait<AssemblerGroup, N>::type;

public:
  template <class... Args> explicit AssemblerGroup(Args &&... args);

  /**
   * @brief Get the Nth member assembler.
   */
  template <std::size_t N> const member_assembler_type<N> &get() const;

  /**
   * @brief Get the Nth member assembler.
   */
  template <std::size_t N> member_assembler_type<N> &get();

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

/**
 * @brief Specialization: The number of assemblers in the group is the number of
 * SingleElementAssemblers-template-parameters.
 */
template <class... SingleElementAssemblers>
struct NumberOfMembersTypeTrait<AssemblerGroup<SingleElementAssemblers...>>
    : std::integral_constant<std::size_t, sizeof...(SingleElementAssemblers)> {
};

/**
 * @brief Specialization: The Nth member is the Nth element of
 * SingleElementAssemblers.
 */
template <std::size_t N, class... SingleElementAssemblers>
struct MemberTypeTrait<AssemblerGroup<SingleElementAssemblers...>, N> {
  using type =
      typename std::tuple_element<N,
                                  std::tuple<SingleElementAssemblers...>>::type;
};
} // namespace assembly
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

namespace ae108 {
namespace assembly {

template <class... SingleElementAssemblers>
template <class... Args>
AssemblerGroup<SingleElementAssemblers...>::AssemblerGroup(Args &&... args)
    : _assemblers(std::forward<Args>(args)...) {}

template <class... SingleElementAssemblers>
template <std::size_t N>
const typename AssemblerGroup<
    SingleElementAssemblers...>::template member_assembler_type<N> &
AssemblerGroup<SingleElementAssemblers...>::get() const {
  return std::get<N>(_assemblers);
}

template <class... SingleElementAssemblers>
template <std::size_t N>
typename AssemblerGroup<
    SingleElementAssemblers...>::template member_assembler_type<N> &
AssemblerGroup<SingleElementAssemblers...>::get() {
  return std::get<N>(_assemblers);
}

} // namespace assembly
} // namespace ae108
