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

#include "ae108/elements/mesh/generate_quadratic_triangle_mesh.h"

namespace ae108 {
namespace elements {
namespace mesh {

Mesh<tensor::Tensor<double, 2>> generate_quadratic_triangle_mesh(
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

    tensor::Tensor<size_type, 2> half_step;
    for (half_step[0] = 0; half_step[0] < 2 * granularity[0] + 1;
         ++half_step[0]) {
      for (half_step[1] = 0; half_step[1] < 2 * granularity[1] + 1;
           ++half_step[1]) {
        positions.push_back({{
            (static_cast<double>(half_step[0]) * size[0]) / 2. /
                static_cast<double>(granularity[0]),
            (static_cast<double>(half_step[1]) * size[1]) / 2. /
                static_cast<double>(granularity[1]),
        }});
      }
    }

    return positions;
  };

  const auto step_to_index = [&](const size_type step_x, const size_type step_y,
                                 const size_type half_step_x = 0,
                                 const size_type half_step_y = 0) {
    return (2 * step_x + half_step_x) * (2 * granularity[1] + 1) + 2 * step_y +
           half_step_y;
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
            step_to_index(step[0], step[1], 1, 0),
            step_to_index(step[0], step[1], 1, 1),
            step_to_index(step[0], step[1], 0, 1),
        }});
        connectivity.push_back({{
            step_to_index(step[0] + 1, step[1]),
            step_to_index(step[0] + 1, step[1] + 1),
            step_to_index(step[0], step[1] + 1),
            step_to_index(step[0] + 1, step[1], 0, 1),
            step_to_index(step[0], step[1] + 1, 1, 0),
            step_to_index(step[0], step[1], 1, 1),
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