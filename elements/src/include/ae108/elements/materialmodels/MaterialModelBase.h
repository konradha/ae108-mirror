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

#include "ae108/elements/tensor/Tensor.h"
#include <limits>

namespace ae108 {
namespace elements {
namespace materialmodels {

template <class SizeType_, class ValueType_, class RealType_,
          SizeType_ Dimension_, SizeType_ DegreesOfFreedom_ = Dimension_>
struct MaterialModelBase {
  using size_type = SizeType_;
  using value_type = ValueType_;
  using real_type = RealType_;

  using Energy = real_type;
  using Time = real_type;

  /**
   * @brief Use this id to specify an unknown displacement gradient id.
   */
  static constexpr size_type unknown_id() noexcept {
    return std::numeric_limits<size_type>::max();
  }

  static constexpr size_type dimension() noexcept { return Dimension_; };

  static constexpr size_type degrees_of_freedom() noexcept {
    return DegreesOfFreedom_;
  };

  /**
   * @brief Displacement gradient $v_{ij}$ in row-major format v[i][j].
   */
  using DisplacementGradient =
      tensor::Tensor<value_type, degrees_of_freedom(), dimension()>;

  /**
   * @brief Stress $P_{ij} = \delta_{v_{ij}} E(v)$ in row-major format P[i][j].
   */
  using Stress = DisplacementGradient;
  using Strain = Stress;

  /**
   * @brief Tangent matrix $C_{ijkl} = \delta_{v_{kl}} \delta_{v_{ij}} E(v)$ in
   * row-major format C[i][j][k][l].
   */
  using TangentMatrix =
      tensor::Tensor<value_type, degrees_of_freedom(), dimension(),
                     degrees_of_freedom(), dimension()>;

protected:
  ~MaterialModelBase() = default;
};
} // namespace materialmodels
} // namespace elements
} // namespace ae108
