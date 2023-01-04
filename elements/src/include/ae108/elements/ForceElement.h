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
#include "ae108/elements/ComputeStiffnessMatrixTrait.h"
#include "ae108/elements/ElementBase.h"
#include "ae108/elements/tensor/as_vector.h"

namespace ae108 {
namespace elements {

/**
 * @brief A single-vertex element that applies a force at that vertex.
 */
template <std::size_t Dimension_, class ValueType_ = double,
          class RealType_ = double>
struct ForceElement final
    : ElementBase<ForceElement<Dimension_, ValueType_, RealType_>, std::size_t,
                  ValueType_, RealType_, 1 /* size */, Dimension_,
                  Dimension_ /* degrees of freedom */> {
public:
  using Force = tensor::Tensor<typename ForceElement::value_type, Dimension_>;

  /**
   * @brief Constructs the element by providing the force to apply.
   */
  explicit ForceElement(Force force) noexcept : force_(std::move(force)) {}

  /**
   * @brief Returns the forces applied at the single vertex.
   */
  const Force &force() const { return force_; }

private:
  Force force_;
};

template <std::size_t Dimension_, class ValueType_, class RealType_>
struct ComputeEnergyTrait<ForceElement<Dimension_, ValueType_, RealType_>> {
  template <class Element>
  typename Element::Energy
  operator()(const Element &element,
             const typename Element::NodalDisplacements &u,
             const typename Element::Time &) const noexcept {
    const auto &force = element.force();
    return tensor::as_vector(&u[0]).dot(tensor::as_vector(&force));
  }
};

template <std::size_t Dimension_, class ValueType_, class RealType_>
struct ComputeForcesTrait<ForceElement<Dimension_, ValueType_, RealType_>> {
  template <class Element>
  typename Element::Forces
  operator()(const Element &element,
             const typename Element::NodalDisplacements &,
             const typename Element::Time &) const noexcept {
    return {element.force()};
  }
};

template <std::size_t Dimension_, class ValueType_, class RealType_>
struct ComputeStiffnessMatrixTrait<
    ForceElement<Dimension_, ValueType_, RealType_>> {
  template <class Element>
  typename Element::StiffnessMatrix
  operator()(const Element &, const typename Element::NodalDisplacements &,
             const typename Element::Time &) const noexcept {
    return Element::StiffnessMatrix::Zero();
  }
};

} // namespace elements
} // namespace ae108