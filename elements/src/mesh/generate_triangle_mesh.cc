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

#include "ae108/elements/mesh/generate_triangle_mesh.h"

namespace ae108 {
namespace elements {
namespace mesh {

Mesh<tensor::Tensor<double, 2>> generate_triangle_mesh(
    const tensor::Tensor<double, 2> &size,
    const tensor::Tensor<std::size_t, 2> &granularity) noexcept {
  using size_type = std::size_t;

  using Point = tensor::Tensor<double, 2>;
  using Mesh = Mesh<Point>;

  if (granularity[0] * granularity[1] == 0) {
    return Mesh({}, {});
  }

  const auto generate_positions = [&]() {
    Mesh::Positions positions;
    positions.reserve((granularity[0] + 1) * (granularity[1] + 1));

    tensor::Tensor<size_type, 2> step;
    for (step[0] = 0; step[0] < granularity[0] + 1; ++step[0]) {
      for (step[1] = 0; step[1] < granularity[1] + 1; ++step[1]) {
        positions.push_back({{
            (static_cast<double>(step[0]) * size[0]) /
                static_cast<double>(granularity[0]),
            (static_cast<double>(step[1]) * size[1]) /
                static_cast<double>(granularity[1]),
        }});
      }
    }

    return positions;
  };

  const auto step_to_index = [&](const size_type step_x,
                                 const size_type step_y) {
    return step_x * (granularity[1] + 1) + step_y;
  };

  const auto generate_connectivity = [&]() {
    Mesh::Connectivity connectivity;
    connectivity.reserve(2 * granularity[0] * granularity[1]);

    tensor::Tensor<size_type, 2> step = {{0, 0}};
    for (step[0] = 0; step[0] < granularity[0]; ++step[0]) {
      for (step[1] = 0; step[1] < granularity[1]; ++step[1]) {
        connectivity.push_back({{
            step_to_index(step[0], step[1]),
            step_to_index(step[0] + 1, step[1]),
            step_to_index(step[0], step[1] + 1),
        }});
        connectivity.push_back({{
            step_to_index(step[0] + 1, step[1]),
            step_to_index(step[0] + 1, step[1] + 1),
            step_to_index(step[0], step[1] + 1),
        }});
      }
    }

    return connectivity;
  };

  return Mesh(generate_connectivity(), generate_positions());
}

} // namespace mesh
} // namespace elements
} // namespace ae108