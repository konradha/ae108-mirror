// © 2021 ETH Zurich, Mechanics and Materials Lab
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
#include <cstddef>
#include <vector>

namespace ae108 {
namespace elements {
namespace mesh {

/**
 * @brief Creates a simple mesh of cuboid shape. First the cuboid is split
 * into smaller cuboids as defined by `granularity`. Each of these smaller
 * cuboids is then split into 5 tetrehedra.
 *
 * @param size The vertex of the cuboid that is opposite to [0., 0., 0.].
 * @param granularity The number of quintets of tetrahedra along the three axes.
 */
Mesh<tensor::Tensor<double, 3>> generate_tetrahedron_mesh(
    const tensor::Tensor<double, 3> &size,
    const tensor::Tensor<std::size_t, 3> &granularity) noexcept;

} // namespace mesh
} // namespace elements
} // namespace ae108