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

#include <Eigen/Dense>
#include <array>
#include <cstddef>
#include <gmock/gmock.h>

namespace ae108 {
namespace assembly {
namespace test {

template <class ValueType_> struct Element {
  static constexpr std::size_t SpatialDimension = 1;
  static constexpr std::size_t DegreesOfFreedom = 1;
  static constexpr std::size_t NumberOfNodes = 1;

  using value_type = ValueType_;

  using Vector = Eigen::Matrix<value_type, 1, 1>;
  using StiffnessMatrix = Vector;
  using Strain = Vector;

  using NodalDisplacements = std::array<Vector, NumberOfNodes>;
  using Forces = NodalDisplacements;

  MOCK_CONST_METHOD2_T(computeEnergy,
                       double(const NodalDisplacements &, const double));
  MOCK_CONST_METHOD2_T(computeForces,
                       Forces(const NodalDisplacements &, const double));
  MOCK_CONST_METHOD2_T(computeStiffnessMatrix,
                       StiffnessMatrix(const NodalDisplacements &,
                                       const double));
  MOCK_CONST_METHOD0_T(computeMassMatrix, StiffnessMatrix());
  MOCK_CONST_METHOD0_T(computeLumpedMassMatrix, StiffnessMatrix());
  MOCK_CONST_METHOD0_T(computeConsistentMassMatrix, StiffnessMatrix());
  MOCK_METHOD2_T(updateInternalVariables,
                 void(const NodalDisplacements &, const double));
};

template <class ValueType_>
auto compute_mass_matrix(const Element<ValueType_> &element) {
  return element.computeMassMatrix();
}

template <class ValueType_>
auto compute_consistent_mass_matrix(const Element<ValueType_> &element) {
  return element.computeConsistentMassMatrix();
}

template <class ValueType_>
auto compute_lumped_mass_matrix(const Element<ValueType_> &element) {
  return element.computeLumpedMassMatrix();
}

} // namespace test
} // namespace assembly
} // namespace ae108
