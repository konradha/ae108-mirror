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

#include "ae108/elements/mesh/generate_tetrahedron_mesh.h"
#include <array>
#include <functional>
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

constexpr auto tetrahedra_per_cuboid = size_type{5};
constexpr auto points_per_tetrahedron = size_type{4};

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
  connectivity.reserve(number_of_cuboids(granularity) * tetrahedra_per_cuboid);

  using ElementVertexSteps =
      std::array<std::array<Index, points_per_tetrahedron>,
                 tetrahedra_per_cuboid>;
  constexpr ElementVertexSteps quintet_A = {{
      {{{{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}}},
      {{{{1, 0, 0}}, {{1, 1, 1}}, {{0, 0, 1}}, {{1, 0, 1}}}},
      {{{{1, 0, 0}}, {{0, 1, 0}}, {{1, 1, 1}}, {{1, 1, 0}}}},
      {{{{1, 1, 1}}, {{0, 1, 0}}, {{0, 0, 1}}, {{0, 1, 1}}}},
      {{{{1, 1, 1}}, {{1, 0, 0}}, {{0, 0, 1}}, {{0, 1, 0}}}},
  }};
  constexpr ElementVertexSteps quintet_B = {{
      {{{{0, 0, 0}}, {{1, 0, 0}}, {{1, 1, 0}}, {{1, 0, 1}}}},
      {{{{1, 1, 0}}, {{0, 1, 1}}, {{1, 0, 1}}, {{1, 1, 1}}}},
      {{{{0, 0, 0}}, {{0, 1, 1}}, {{1, 1, 0}}, {{0, 1, 0}}}},
      {{{{0, 0, 0}}, {{1, 0, 1}}, {{0, 1, 1}}, {{0, 0, 1}}}},
      {{{{1, 0, 1}}, {{0, 1, 1}}, {{1, 1, 0}}, {{0, 0, 0}}}},
  }};

  const auto add_quintet = [&](const ElementVertexSteps &quintet,
                               const Index &offset) {
    for (auto &&element : quintet) {
      connectivity.emplace_back();
      connectivity.back().reserve(points_per_tetrahedron);
      for (auto &&steps : element) {
        connectivity.back().push_back(steps_to_id({{
                                                      steps[0] + offset[0],
                                                      steps[1] + offset[1],
                                                      steps[2] + offset[2],
                                                  }},
                                                  granularity));
      }
    }
  };

  auto offset = Index();
  for (offset[2] = 0; offset[2] < granularity[2]; ++offset[2])
    for (offset[1] = 0; offset[1] < granularity[1]; ++offset[1])
      for (offset[0] = 0; offset[0] < granularity[0]; ++offset[0]) {
        const auto &quintet =
            (offset[0] + offset[1] + offset[2]) % 2 ? quintet_B : quintet_A;
        add_quintet(quintet, offset);
      }

  return connectivity;
}

} // namespace

Mesh<Point> generate_tetrahedron_mesh(const Point &size,
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