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

#include "ae108/elements/AutomaticForcesTrait.h"
#include "ae108/elements/AutomaticStiffnessMatrixTrait.h"
#include "ae108/elements/ComputeEnergyTrait.h"
#include "ae108/elements/ComputeForcesTrait.h"
#include "ae108/elements/ComputeStiffnessMatrixTrait.h"
#include "ae108/elements/ElementBase.h"

namespace ae108 {
namespace elements {

/**
 * @brief A minimal example of an element. Note that only the energy is
 * defined.
 */
template <std::size_t Size_, std::size_t Dimension_,
          std::size_t DegreesOfFreedom_ = Dimension_>
struct Minimal final
    : ElementBase<Minimal<Size_, Dimension_, DegreesOfFreedom_>, std::size_t,
                  double, double, Size_, Dimension_, DegreesOfFreedom_> {};

/**
 * @brief Always returns 0.
 */
template <std::size_t Size_, std::size_t Dimension_,
          std::size_t DegreesOfFreedom_>
struct ComputeEnergyTrait<Minimal<Size_, Dimension_, DegreesOfFreedom_>> {
  template <class Element>
  typename Element::Energy
  operator()(const Element &, const typename Element::NodalDisplacements &,
             const typename Element::Time &) const noexcept {
    return 0.;
  }
};

/**
 * @brief Computes the forces by differentiating the energy.
 */
template <std::size_t Size_, std::size_t Dimension_,
          std::size_t DegreesOfFreedom_>
struct ComputeForcesTrait<Minimal<Size_, Dimension_, DegreesOfFreedom_>>
    : AutomaticForcesTrait<Minimal<Size_, Dimension_, DegreesOfFreedom_>> {};

/**
 * @brief Computes the stiffness matrix by differentiating the forces.
 */
template <std::size_t Size_, std::size_t Dimension_,
          std::size_t DegreesOfFreedom_>
struct ComputeStiffnessMatrixTrait<
    Minimal<Size_, Dimension_, DegreesOfFreedom_>>
    : AutomaticStiffnessMatrixTrait<
          Minimal<Size_, Dimension_, DegreesOfFreedom_>> {};

} // namespace elements
} // namespace ae108