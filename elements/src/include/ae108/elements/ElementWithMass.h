// © 2020 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/elements/ComputeEnergyTrait.h"
#include "ae108/elements/ComputeForcesTrait.h"
#include "ae108/elements/ComputeMassMatrixTrait.h"
#include "ae108/elements/ComputeStiffnessMatrixTrait.h"
#include "ae108/elements/ElementBase.h"
#include "ae108/elements/compute_energy.h"
#include "ae108/elements/compute_forces.h"
#include "ae108/elements/compute_stiffness_matrix.h"
#include <utility>

namespace ae108 {
namespace elements {

/**
 * @brief Extends a given element to also store a mass matrix.
 */
template <class Element_>
struct ElementWithMass final
    : ElementBase<ElementWithMass<Element_>, typename Element_::size_type,
                  typename Element_::value_type, typename Element_::real_type,
                  Element_::size(), Element_::dimension(),
                  Element_::degrees_of_freedom()> {
public:
  using Element = Element_;
  using MassMatrix = typename ElementWithMass::StiffnessMatrix;

  explicit ElementWithMass(Element element, MassMatrix mass_matrix) noexcept
      : element_(std::move(element)), mass_matrix_(std::move(mass_matrix)) {}

  const MassMatrix &mass_matrix() const & { return mass_matrix_; }
  const Element &unwrap() const & { return element_; }

private:
  Element element_;
  MassMatrix mass_matrix_;
};

template <class Element_> struct ComputeEnergyTrait<ElementWithMass<Element_>> {
  template <class Element, class... Args>
  auto operator()(const Element &element, Args &&...args) const noexcept {
    return compute_energy(element.unwrap(), std::forward<Args>(args)...);
  }
};

template <class Element_> struct ComputeForcesTrait<ElementWithMass<Element_>> {
  template <class Element, class... Args>
  auto operator()(const Element &element, Args &&...args) const noexcept {
    return compute_forces(element.unwrap(), std::forward<Args>(args)...);
  }
};

template <class Element_>
struct ComputeStiffnessMatrixTrait<ElementWithMass<Element_>> {
  template <class Element, class... Args>
  auto operator()(const Element &element, Args &&...args) const noexcept {
    return compute_stiffness_matrix(element.unwrap(),
                                    std::forward<Args>(args)...);
  }
};

template <class Element_>
struct ComputeMassMatrixTrait<ElementWithMass<Element_>> {
  template <class Element>
  auto &operator()(const Element &element) const &noexcept {
    return element.mass_matrix();
  }
};
} // namespace elements
} // namespace ae108