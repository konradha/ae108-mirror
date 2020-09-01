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

#include <Eigen/Dense>
#include <array>
#include <cstddef>
#include <gmock/gmock.h>

namespace ae108 {
namespace assembly {
namespace test {

struct Element {
  static constexpr std::size_t SpatialDimension = 1;
  static constexpr std::size_t DegreesOfFreedom = 1;
  static constexpr std::size_t NumberOfNodes = 1;

  using Vector = Eigen::Matrix<double, 1, 1>;
  using StiffnessMatrix = Vector;
  using Strain = Vector;

  using NodalDisplacements = std::array<Vector, NumberOfNodes>;
  using Forces = NodalDisplacements;

  MOCK_CONST_METHOD2(computeEnergy,
                     double(const NodalDisplacements &, const double));
  MOCK_CONST_METHOD2(computeForces,
                     Forces(const NodalDisplacements &, const double));
  MOCK_CONST_METHOD2(computeStiffnessMatrix,
                     StiffnessMatrix(const NodalDisplacements &, const double));
  MOCK_CONST_METHOD0(computeLumpedMassMatrix, StiffnessMatrix());
  MOCK_METHOD2(updateInternalVariables,
               void(const NodalDisplacements &, const double));
};
} // namespace test
} // namespace assembly
} // namespace ae108
