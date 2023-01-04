// © 2020 ETH Zurich, Mechanics and Materials Lab
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