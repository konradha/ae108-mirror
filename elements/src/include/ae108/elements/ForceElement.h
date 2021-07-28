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
#include "ae108/elements/ComputeStiffnessMatrixTrait.h"
#include "ae108/elements/ElementBase.h"
#include "ae108/elements/tensor/as_vector.h"

namespace ae108 {
namespace elements {

/**
 * @brief A single-vertex element that applies a force at that vertex.
 */
template <std::size_t DegreesOfFreedom_>
struct ForceElement final
    : ElementBase<ForceElement<DegreesOfFreedom_>, std::size_t, double,
                  1 /* size */, DegreesOfFreedom_> {
public:
  using Force =
      tensor::Tensor<typename ForceElement::value_type, DegreesOfFreedom_>;

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

template <std::size_t DegreesOfFreedom_>
struct ComputeEnergyTrait<ForceElement<DegreesOfFreedom_>> {
  template <class Element>
  typename Element::Energy
  operator()(const Element &element,
             const typename Element::NodalDisplacements &u,
             const typename Element::Time &) const noexcept {
    const auto &force = element.force();
    return tensor::as_vector(&u[0]).dot(tensor::as_vector(&force));
  }
};

template <std::size_t DegreesOfFreedom_>
struct ComputeForcesTrait<ForceElement<DegreesOfFreedom_>> {
  template <class Element>
  typename Element::Forces
  operator()(const Element &element,
             const typename Element::NodalDisplacements &,
             const typename Element::Time &) const noexcept {
    return {element.force()};
  }
};

template <std::size_t DegreesOfFreedom_>
struct ComputeStiffnessMatrixTrait<ForceElement<DegreesOfFreedom_>> {
  template <class Element>
  typename Element::StiffnessMatrix
  operator()(const Element &, const typename Element::NodalDisplacements &,
             const typename Element::Time &) const noexcept {
    return Element::StiffnessMatrix::Zero();
  }
};

} // namespace elements
} // namespace ae108