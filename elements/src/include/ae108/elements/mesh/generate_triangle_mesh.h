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

#include "ae108/elements/mesh/Mesh.h"
#include "ae108/elements/tensor/Tensor.h"
#include <vector>

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
 * @param size The top-right coordinate of the rectangle. The bottom-left
 * coordinate of the rectangle is [0., 0.].
 * @param granularity The number of pairs of triangles along the two axes.
 */
Mesh<tensor::Tensor<double, 2>> generate_triangle_mesh(
    const tensor::Tensor<double, 2> &size,
    const tensor::Tensor<std::size_t, 2> &granularity) noexcept;

} // namespace mesh
} // namespace elements
} // namespace ae108