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