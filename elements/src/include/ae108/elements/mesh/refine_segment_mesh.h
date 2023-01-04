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
                                const double max_segment_length) noexcept;

extern template Mesh<tensor::Tensor<double, 2>>
refine_segment_mesh(const Mesh<tensor::Tensor<double, 2>> &unrefined_mesh,
                    const double max_segment_length) noexcept;

extern template Mesh<tensor::Tensor<double, 3>>
refine_segment_mesh(const Mesh<tensor::Tensor<double, 3>> &unrefined_mesh,
                    const double max_segment_length) noexcept;

} // namespace mesh
} // namespace elements
} // namespace ae108

namespace ae108 {
namespace elements {
namespace mesh {

template <class Point>
Mesh<Point> refine_segment_mesh(const Mesh<Point> &unrefined_mesh,
                                const double max_segment_length) noexcept {
  using Mesh = Mesh<Point>;
  typename Mesh::Positions positions;
  typename Mesh::Connectivity connectivity;

  positions.reserve(unrefined_mesh.number_of_positions());
  connectivity.reserve(unrefined_mesh.connectivity().size());

  const auto number_of_positions = unrefined_mesh.number_of_positions();
  for (auto i = decltype(number_of_positions){0}; i < number_of_positions; i++)
    positions.push_back(unrefined_mesh.position_of_vertex(i));

  for (const auto &segment : unrefined_mesh.connectivity()) {
    assert(segment.size() == 2u && "A segment connects two vertices.");

    const auto &source =
        tensor::as_vector(&unrefined_mesh.position_of_vertex(segment[0]));

    const auto &target =
        tensor::as_vector(&unrefined_mesh.position_of_vertex(segment[1]));

    const auto divisions =
        std::ceil((target - source).norm() / max_segment_length);

    auto source_id = segment[0];
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
