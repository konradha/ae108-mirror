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

#include "ae108/elements/compute_energy.h"
#include "ae108/elements/compute_forces.h"
#include "ae108/elements/compute_stiffness_matrix.h"
#include "ae108/elements/tensor/Tensor.h"
#include <Eigen/Core>
#include <cstddef>

namespace ae108 {
namespace elements {

template <class Derived_, class SizeType_, class ValueType_, SizeType_ Size_,
          SizeType_ DegreesOfFreedom_>
struct ElementBase {
  using size_type = SizeType_;
  using value_type = ValueType_;

  using Energy = value_type;
  using Time = value_type;

  /**
   * @brief The displacements per node.
   */
  using NodalDisplacements =
      tensor::Tensor<value_type, Size_, DegreesOfFreedom_>;

  /**
   * @brief The forces equal to $d_{ij} E$.
   * Here, d_ij refers to the derivative with respect to jth degree of freedom
   * of the ith node.
   */
  using Forces = NodalDisplacements;

  /**
   * @brief The stiffness matrix equal to $d_{ij} d_{kl} E$.
   * Here, d_ij refers to the derivative with respect to jth degree of freedom
   * of the ith node.
   */
  using StiffnessMatrix =
      Eigen::Matrix<value_type, Size_ * DegreesOfFreedom_,
                    Size_ * DegreesOfFreedom_, Eigen::RowMajor>;

  /**
   * @brief Number of element nodes / shape functions.
   */
  static constexpr size_type size() noexcept { return Size_; }

  /**
   * @brief Number of degrees of freedom.
   */
  static constexpr size_type degrees_of_freedom() noexcept {
    return DegreesOfFreedom_;
  }

  /**
   * @brief Computes the energy for the given displacements.
   * @remark Delegates to the corresponding free function.
   */
  Energy computeEnergy(const NodalDisplacements &displacements,
                       const Time &time) const {
    return compute_energy(static_cast<const Derived_ &>(*this), displacements,
                          time);
  }

  /**
   * @brief Computes the forces for the given displacements.
   * @remark Delegates to the corresponding free function.
   */
  Forces computeForces(const NodalDisplacements &displacements,
                       const Time &time) const {
    return compute_forces(static_cast<const Derived_ &>(*this), displacements,
                          time);
  }

  /**
   * @brief Computes the stiffness matrix for the given displacements.
   * @remark Delegates to the corresponding free function.
   */
  StiffnessMatrix
  computeStiffnessMatrix(const NodalDisplacements &displacements,
                         const Time &time) const {
    return compute_stiffness_matrix(static_cast<const Derived_ &>(*this),
                                    displacements, time);
  }

protected:
  ~ElementBase() = default;
};

} // namespace elements
} // namespace ae108
