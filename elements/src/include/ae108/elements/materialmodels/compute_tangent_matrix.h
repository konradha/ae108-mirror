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

#include "ae108/elements/materialmodels/ComputeTangentMatrixTrait.h"

namespace ae108 {
namespace elements {
namespace materialmodels {

template <class MaterialModel_>
typename MaterialModel_::TangentMatrix compute_tangent_matrix(
    const MaterialModel_ &model, const typename MaterialModel_::size_type id,
    const typename MaterialModel_::DisplacementGradient &gradient,
    const typename MaterialModel_::Time time) noexcept {
  return ComputeTangentMatrixTrait<MaterialModel_>()(model, id, gradient, time);
}

} // namespace materialmodels
} // namespace elements
} // namespace ae108
