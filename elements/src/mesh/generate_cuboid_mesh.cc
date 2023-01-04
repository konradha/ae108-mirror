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

#include "ae108/elements/mesh/generate_cuboid_mesh.h"
#include <array>
#include <numeric>

namespace ae108 {
namespace elements {
namespace mesh {

namespace {

using Point = Mesh<tensor::Tensor<double, 3>>::Point;
using size_type = Mesh<Point>::size_type;
using value_type = Point::value_type;
using Positions = Mesh<Point>::Positions;
using Connectivity = Mesh<Point>::Connectivity;
using Index = tensor::Tensor<size_type, 3>;

constexpr auto points_per_cuboid = size_type{8};

/**
 * @brief Computes the number of smaller cuboids the cuboid
 * will be split into.
 */
size_type number_of_cuboids(const Index &granularity) noexcept {
  return std::accumulate(granularity.begin(), granularity.end(), size_type{1},
                         std::multiplies<size_type>{});
}

/**
 * @brief Converts the steps in x, y, and z
 * direction to the ID of the point at that location.
 */
size_type steps_to_id(const Index &steps, const Index &granularity) noexcept {
  return steps[2] * (granularity[1] + size_type{1}) *
             (granularity[0] + size_type{1}) +
         steps[1] * (granularity[0] + size_type{1}) + steps[0];
}

/**
 * @brief Converts the steps in x, y, and z
 * direction to the position of the point at that location.
 */
Point steps_to_position(const Index &steps, const Point &size,
                        const Index &granularity) noexcept {
  return {{
      static_cast<value_type>(steps[0]) * static_cast<value_type>(size[0]) /
          static_cast<value_type>(granularity[0]),
      static_cast<value_type>(steps[1]) * static_cast<value_type>(size[1]) /
          static_cast<value_type>(granularity[1]),
      static_cast<value_type>(steps[2]) * static_cast<value_type>(size[2]) /
          static_cast<value_type>(granularity[2]),
  }};
}

/**
 * @brief Generates the positions of the points in the mesh.
 */
Positions generate_positions(const Point &size,
                             const Index &granularity) noexcept {
  auto positions = Positions();
  positions.reserve(number_of_cuboids(granularity));

  auto steps = Index();
  for (steps[2] = 0; steps[2] <= granularity[2]; ++steps[2])
    for (steps[1] = 0; steps[1] <= granularity[1]; ++steps[1])
      for (steps[0] = 0; steps[0] <= granularity[0]; ++steps[0])
        positions.emplace_back(steps_to_position(steps, size, granularity));

  return positions;
}

/**
 * @brief Generates the connectivity of the mesh.
 */
Connectivity generate_connectivity(const Index &granularity) noexcept {
  auto connectivity = Connectivity();
  connectivity.reserve(number_of_cuboids(granularity));

  constexpr std::array<Index, points_per_cuboid> steps = {{
      {{0, 0, 0}},
      {{1, 0, 0}},
      {{1, 1, 0}},
      {{0, 1, 0}},
      {{0, 0, 1}},
      {{1, 0, 1}},
      {{1, 1, 1}},
      {{0, 1, 1}},
  }};

  const auto add_cuboid = [&](const Index &offset) {
    connectivity.emplace_back();
    connectivity.back().reserve(points_per_cuboid);
    for (auto &&step : steps) {
      connectivity.back().push_back(steps_to_id({{
                                                    step[0] + offset[0],
                                                    step[1] + offset[1],
                                                    step[2] + offset[2],
                                                }},
                                                granularity));
    }
  };

  auto offset = Index();
  for (offset[2] = 0; offset[2] < granularity[2]; ++offset[2])
    for (offset[1] = 0; offset[1] < granularity[1]; ++offset[1])
      for (offset[0] = 0; offset[0] < granularity[0]; ++offset[0])
        add_cuboid(offset);

  return connectivity;
}

} // namespace

Mesh<Point> generate_cuboid_mesh(const Point &size,
                                 const Index &granularity) noexcept {
  if (number_of_cuboids(granularity) == 0) {
    return Mesh<Point>({}, {});
  }

  return Mesh<Point>(generate_connectivity(granularity),
                     generate_positions(size, granularity));
}

} // namespace mesh
} // namespace elements
} // namespace ae108