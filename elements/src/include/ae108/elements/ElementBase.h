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

#include "ae108/elements/compute_energy.h"
#include "ae108/elements/compute_forces.h"
#include "ae108/elements/compute_stiffness_matrix.h"
#include "ae108/elements/tensor/Tensor.h"
#include <Eigen/Core>
#include <cstddef>

namespace ae108 {
namespace elements {

template <class Derived_, class SizeType_, class ValueType_, class RealType_,
          SizeType_ Size_, SizeType_ Dimension_, SizeType_ DegreesOfFreedom_>
struct ElementBase {
  using size_type = SizeType_;
  using value_type = ValueType_;
  using real_type = RealType_;

  using Energy = real_type;
  using Time = real_type;

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
   * @brief The dimension of physical space.
   */
  static constexpr size_type dimension() noexcept { return Dimension_; }

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
