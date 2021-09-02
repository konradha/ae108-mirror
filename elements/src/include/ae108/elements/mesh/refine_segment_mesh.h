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
#include "ae108/elements/tensor/as_vector.h"

namespace ae108 {
namespace elements {
namespace mesh {

/**
 * @brief Refines a given mesh of segments. Each segment is refined into smaller
 * segments of equal length such that no segment is longer than
 * `max_segment_length`.
 *
 * @param unrefined_mesh The unrefined segment mesh.
 * @param max_segment_length The maximum allowed segment length.
 */
template <class Point>
Mesh<Point> refine_segment_mesh(const Mesh<Point> &unrefined_mesh,
                                const double &max_segment_length) noexcept;

extern template Mesh<tensor::Tensor<double, 2>>
refine_segment_mesh(const Mesh<tensor::Tensor<double, 2>> &unrefined_mesh,
                    const double &max_segment_length) noexcept;

extern template Mesh<tensor::Tensor<double, 3>>
refine_segment_mesh(const Mesh<tensor::Tensor<double, 3>> &unrefined_mesh,
                    const double &max_segment_length) noexcept;

} // namespace mesh
} // namespace elements
} // namespace ae108

namespace ae108 {
namespace elements {
namespace mesh {

template <class Point>
Mesh<Point> refine_segment_mesh(const Mesh<Point> &unrefined_mesh,
                                const double &max_segment_length) noexcept {
  using Mesh = Mesh<Point>;
  typename Mesh::Positions positions;
  typename Mesh::Connectivity connectivity;

  for (std::size_t i = 0; i < unrefined_mesh.number_of_positions(); i++)
    positions.push_back(unrefined_mesh.position_of_vertex(i));

  for (const auto &segment : unrefined_mesh.connectivity()) {

    const auto &source =
        tensor::as_vector(&unrefined_mesh.position_of_vertex(segment[0]));

    const auto &target =
        tensor::as_vector(&unrefined_mesh.position_of_vertex(segment[1]));

    const auto divisions =
        std::ceil((target - source).norm() / max_segment_length);

    std::size_t source_id = segment[0];
    for (std::size_t increment = 1; increment < divisions; increment++) {
      connectivity.push_back({source_id, positions.size()});
      positions.push_back([&]() {
        Point point;
        tensor::as_vector(&point) =
            source + increment * (target - source) / divisions;
        return point;
      }());
      source_id = positions.size() - 1;
    }
    connectivity.push_back({source_id, segment[1]});
  }
  return Mesh(connectivity, positions);
}

} // namespace mesh
} // namespace elements
} // namespace ae108
