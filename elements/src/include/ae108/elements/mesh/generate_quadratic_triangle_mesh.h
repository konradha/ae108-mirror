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

#include "ae108/elements/mesh/Mesh.h"
#include "ae108/elements/tensor/Tensor.h"

namespace ae108 {
namespace elements {
namespace mesh {

/**
 * @brief Creates a simple triangulation of a rectangle. For instance, it
 * creates the following triangulation for granularity [2, 1]:
 *
 *  *-----*-----*
 *  *  \  |  \  |
 *  *-----*-----*
 *
 * The function generates quadratic triangles. In particular, it uses
 * the following order for the nodes:
 *
 * 2
 * | \
 * 5   4
 * |    \
 * 0--3--1
 *
 * @param size The top-right coordinate of the rectangle. The bottom-left
 * coordinate of the rectangle is [0., 0.].
 * @param granularity The number of pairs of triangles along the two axes.
 */
Mesh<tensor::Tensor<double, 2>> generate_quadratic_triangle_mesh(
    const tensor::Tensor<double, 2> &size,
    const tensor::Tensor<std::size_t, 2> &granularity) noexcept;

} // namespace mesh
} // namespace elements
} // namespace ae108