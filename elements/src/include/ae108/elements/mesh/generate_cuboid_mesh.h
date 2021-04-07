// © 2021 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/elements/mesh/Mesh.h"
#include "ae108/elements/tensor/Tensor.h"
#include <cstddef>
#include <vector>

namespace ae108 {
namespace elements {
namespace mesh {

/**
 * @brief Creates a simple mesh of cuboid shape by splitting
 * into smaller cuboids as defined by `granularity`.
 *
 * @param size The vertex of the cuboid that is opposite to [0., 0., 0.].
 * @param granularity The number of cuboids along the three axes.
 */
Mesh<tensor::Tensor<double, 3>> generate_cuboid_mesh(
    const tensor::Tensor<double, 3> &size,
    const tensor::Tensor<std::size_t, 3> &granularity) noexcept;

} // namespace mesh
} // namespace elements
} // namespace ae108