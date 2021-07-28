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
#include "ae108/assembly/DefaultFeaturePlugins.h"
#include "ae108/assembly/FeaturePlugin.h"
#include "ae108/cpppetsc/LocalElementView.h"
#include "ae108/cpppetsc/Mesh_fwd.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"
#include <deque>
#include <range/v3/view/subrange.hpp>
#include <utility>

namespace ae108 {
namespace assembly {

/**
 * @brief A cpppetsc-Mesh-based assembler.
 *
 * @tparam Element The element type.
 * @tparam Plugins A FeaturePlugins<...> list of feature plugins.
 */
template <class Element, class Plugins = DefaultFeaturePlugins,
          class Policy = cpppetsc::SequentialComputePolicy>
class Assembler
    : public Plugins::template type<Assembler<Element, Plugins, Policy>> {
  using PluginBase =
      typename Plugins::template type<Assembler<Element, Plugins, Policy>>;

  template <class Assembler, template <class> class Plugin>
  friend class FeaturePlugin;

  /**
   * @brief A class that is essentially a pair of an ElementView and an Element
   * instance.
   */
  class AnnotatedElement;

  using ElementContainer = std::deque<AnnotatedElement>;
  using iterator = typename ElementContainer::iterator;
  using const_iterator = typename ElementContainer::const_iterator;

public:
  using mesh_type = typename MeshTypeTrait<Assembler>::type;
  using element_type = typename ElementTypeTrait<Assembler>::type;

  using size_type = typename SizeTypeTrait<Assembler>::type;
  using vector_type = typename VectorTypeTrait<Assembler>::type;
  using matrix_type = typename MatrixTypeTrait<Assembler>::type;
  using value_type = typename ValueTypeTrait<Assembler>::type;

  /**
   * @param args Are used to construct the plugins.
   */
  template <class... Args> explicit Assembler(Args &&... args);

  using ElementView = typename cpppetsc::LocalElementView<mesh_type>;

  /**
   * @brief Emplace (i.e. construct directly into the internal container without
   * moving) a local Element instance.
   *
   * @remark Use methods of the mesh to know which elements are local, and their
   * vertices.
   *
   * @param view A view of the mesh element. A copy of this view will be stored.
   * @param constructorArguments The parameters that are passed to the Element
   * constructor.
   */
  template <class... Args>
  void emplaceElement(ElementView view, Args &&... constructorArguments);

  /**
   * @brief Returns a range of iterators pointing to a struct with two methods:
   * instance() and meshView().
   */
  auto meshElements();

  /**
   * @brief Returns a range of iterators pointing to a struct with two methods:
   * instance() and meshView().
   */
  auto meshElements() const;

private:
  ElementContainer _elements;
};

/**
 * @brief Deduces the Policy parameter of Assembler.
 */
template <class Element, class Plugins, class Policy>
struct PolicyTypeTrait<Assembler<Element, Plugins, Policy>> {
  using type = Policy;
};

/**
 * @brief Type is an alias for the mesh type used by Assembler.
 */
template <class Element, class Plugins, class Policy>
struct MeshTypeTrait<Assembler<Element, Plugins, Policy>> {
  using type = typename cpppetsc::Mesh<
      typename PolicyTypeTrait<Assembler<Element, Plugins, Policy>>::type>;
};

/**
 * @brief Type is an alias for the vector type used by Assembler.
 */
template <class Element, class Plugins, class Policy>
struct VectorTypeTrait<Assembler<Element, Plugins, Policy>> {
  using type = typename MeshTypeTrait<
      Assembler<Element, Plugins, Policy>>::type::vector_type;
};

/**
 * @brief Type is an alias for the matrix type used by Assembler.
 */
template <class Element, class Plugins, class Policy>
struct MatrixTypeTrait<Assembler<Element, Plugins, Policy>> {
  using type = typename MeshTypeTrait<
      Assembler<Element, Plugins, Policy>>::type::matrix_type;
};

/**
 * @brief Type is an alias for the size type used by Assembler.
 */
template <class Element, class Plugins, class Policy>
struct SizeTypeTrait<Assembler<Element, Plugins, Policy>> {
  using type = typename MeshTypeTrait<
      Assembler<Element, Plugins, Policy>>::type::size_type;
};

/**
 * @brief Type is an alias for the value type used by Assembler.
 */
template <class Element, class Plugins, class Policy>
struct ValueTypeTrait<Assembler<Element, Plugins, Policy>> {
  using type = typename MeshTypeTrait<
      Assembler<Element, Plugins, Policy>>::type::value_type;
};

/**
 * @brief Type is an alias for the element type used by Assembler.
 */
template <class Element, class Plugins, class Policy>
struct ElementTypeTrait<Assembler<Element, Plugins, Policy>> {
  using type = Element;
};

/**
 * @brief Type is an alias for the plugin type used by Assembler.
 */
template <class Element, class Plugins, class Policy>
struct PluginTypeTrait<Assembler<Element, Plugins, Policy>> {
  using type = Plugins;
};
} // namespace assembly
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

#include "ae108/cpppetsc/Mesh.h"

namespace ae108 {
namespace assembly {

template <class Element, class Plugins, class Policy>
template <class... Args>
Assembler<Element, Plugins, Policy>::Assembler(Args &&... args)
    : PluginBase(std::forward<Args>(args)...) {}

template <class Element, class Plugins, class Policy>
template <class... Args>
void Assembler<Element, Plugins, Policy>::emplaceElement(
    ElementView view, Args &&... constructorArguments) {
  _elements.emplace_back(std::move(view),
                         std::forward<Args>(constructorArguments)...);
}

template <class Element, class Plugins, class Policy>
auto Assembler<Element, Plugins, Policy>::meshElements() {
  return ranges::make_subrange(_elements.begin(), _elements.end());
}

template <class Element, class Plugins, class Policy>
auto Assembler<Element, Plugins, Policy>::meshElements() const {
  return ranges::make_subrange(_elements.begin(), _elements.end());
}

template <class Element, class Plugins, class Policy>
class Assembler<Element, Plugins, Policy>::AnnotatedElement {
public:
  /**
   * @param view A copy of this view is stored.
   * @param args Are used to construct in-place the element instance.
   */
  template <class... Args>
  explicit AnnotatedElement(ElementView view, Args &&... args)
      : _meshView(std::move(view)), _instance(std::forward<Args>(args)...) {}

  const ElementView &meshView() const { return _meshView; }

  const Element &instance() const { return _instance; }
  Element &instance() { return _instance; }

private:
  ElementView _meshView;
  Element _instance;
};
} // namespace assembly
} // namespace ae108
