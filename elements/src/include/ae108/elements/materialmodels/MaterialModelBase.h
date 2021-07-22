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

#include "ae108/elements/tensor/Tensor.h"
#include <limits>

namespace ae108 {
namespace elements {
namespace materialmodels {

template <class SizeType_, class ValueType_, SizeType_ Dimension_,
          SizeType_ DegreesOfFreedom_ = Dimension_>
struct MaterialModelBase {
  using size_type = SizeType_;
  using value_type = ValueType_;

  using Energy = value_type;
  using Time = value_type;

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
