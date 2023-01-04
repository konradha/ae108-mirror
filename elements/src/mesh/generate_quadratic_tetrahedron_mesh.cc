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

#include "ae108/elements/mesh/generate_quadratic_tetrahedron_mesh.h"
#include <array>
#include <cassert>
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
constexpr auto points_per_quadratic_tetrahedron = size_type{10};

/**
 * @brief Scales the provided `index` by a factor of 2.
 */
constexpr Index scale(const Index &index) noexcept {
  return {{
      size_type{2} * index[0],
      size_type{2} * index[1],
      size_type{2} * index[2],
  }};
}

/**
 * @brief Adds the provided indices.
 */
constexpr Index add(const Index &lhs, const Index &rhs) noexcept {
  return {{
      lhs[0] + rhs[0],
      lhs[1] + rhs[1],
      lhs[2] + rhs[2],
  }};
}

/**
 * @brief Returns true if and only if `half_steps` points to the
 * center of a cuboid.
 */
constexpr bool is_center(const Index &half_steps) noexcept {
  return half_steps[0] % size_type{2} && half_steps[1] % size_type{2} &&
         half_steps[2] % size_type{2};
}

/**
 * @brief Computes the number of smaller cuboids the cuboid
 * will be split into.
 */
size_type number_of_cuboids(const Index &granularity) noexcept {
  return std::accumulate(granularity.begin(), granularity.end(), size_type{1},
                         std::multiplies<size_type>{});
}

/**
 * @brief Converts the half steps in x, y, and z
 * direction to the ID of the point at that location.
 *
 * @pre `half_steps` does not point to the center of a cuboid.
 */
size_type half_steps_to_id(const Index &half_steps,
                           const Index &granularity) noexcept {
  assert(!is_center(half_steps) &&
         "The point in the center of the cuboid is not used.");

  const auto scaled_granularity = scale(granularity);

  const auto x_correction =
      (half_steps[1] % size_type{2} && half_steps[2] % size_type{2})
          ? half_steps[0] / size_type{2}
          : size_type{0};
  const auto x_contribution = half_steps[0] - x_correction;

  const auto y_correction =
      (half_steps[2] % 2) ? (half_steps[1] / size_type{2}) * granularity[0]
                          : size_type{0};
  const auto y_contribution =
      half_steps[1] * (scaled_granularity[0] + size_type{1}) - y_correction;

  const auto z_correction =
      (half_steps[2] / size_type{2}) * granularity[0] * granularity[1];
  const auto z_contribution = half_steps[2] *
                                  (scaled_granularity[1] + size_type{1}) *
                                  (scaled_granularity[0] + size_type{1}) -
                              z_correction;

  return x_contribution + y_contribution + z_contribution;
}

/**
 * @brief Converts the half steps in x, y, and z
 * direction to the position of the point at that location.
 */
Point half_steps_to_position(const Index &half_steps, const Point &size,
                             const Index &granularity) {
  const auto scaled_granularity = scale(granularity);
  return {{
      static_cast<value_type>(half_steps[0]) *
          static_cast<value_type>(size[0]) /
          static_cast<value_type>(scaled_granularity[0]),
      static_cast<value_type>(half_steps[1]) *
          static_cast<value_type>(size[1]) /
          static_cast<value_type>(scaled_granularity[1]),
      static_cast<value_type>(half_steps[2]) *
          static_cast<value_type>(size[2]) /
          static_cast<value_type>(scaled_granularity[2]),
  }};
}

/**
 * @brief Generates the positions of the points in the mesh.
 */
Positions generate_positions(const Point &size,
                             const Index &granularity) noexcept {
  const auto scaled_granularity = scale(granularity);

  auto positions = Positions();
  positions.reserve(number_of_cuboids(granularity));

  auto half_steps = Index();
  for (half_steps[2] = 0; half_steps[2] <= scaled_granularity[2];
       ++half_steps[2])
    for (half_steps[1] = 0; half_steps[1] <= scaled_granularity[1];
         ++half_steps[1])
      for (half_steps[0] = 0; half_steps[0] <= scaled_granularity[0];
           ++half_steps[0]) {
        if (!is_center(half_steps)) {
          positions.emplace_back(
              half_steps_to_position(half_steps, size, granularity));
        }
      }

  return positions;
}

template <size_type Size_>
using ElementVertexSteps =
    std::array<std::array<Index, Size_>, tetrahedra_per_cuboid>;

/**
 * @brief Converts the description of a linear tetrahedron to the description
 * of a quadratic tetrahedron.
 */
constexpr ElementVertexSteps<points_per_quadratic_tetrahedron>::value_type
to_quadratic_tetrahedron(
    const ElementVertexSteps<points_per_tetrahedron>::value_type
        &steps) noexcept {
  return {{
      scale(steps[0]),
      scale(steps[1]),
      scale(steps[2]),
      scale(steps[3]),
      add(steps[0], steps[1]),
      add(steps[1], steps[2]),
      add(steps[2], steps[0]),
      add(steps[0], steps[3]),
      add(steps[1], steps[3]),
      add(steps[2], steps[3]),
  }};
}

/**
 * @brief Generates the connectivity of the mesh.
 */
Connectivity generate_connectivity(const Index &granularity) noexcept {
  auto connectivity = Connectivity();
  connectivity.reserve(number_of_cuboids(granularity) * tetrahedra_per_cuboid);

  constexpr ElementVertexSteps<points_per_quadratic_tetrahedron> quintet_A = {{
      to_quadratic_tetrahedron(
          {{{{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}}}),
      to_quadratic_tetrahedron(
          {{{{1, 0, 0}}, {{1, 1, 1}}, {{0, 0, 1}}, {{1, 0, 1}}}}),
      to_quadratic_tetrahedron(
          {{{{1, 0, 0}}, {{0, 1, 0}}, {{1, 1, 1}}, {{1, 1, 0}}}}),
      to_quadratic_tetrahedron(
          {{{{1, 1, 1}}, {{0, 1, 0}}, {{0, 0, 1}}, {{0, 1, 1}}}}),
      to_quadratic_tetrahedron(
          {{{{1, 1, 1}}, {{1, 0, 0}}, {{0, 0, 1}}, {{0, 1, 0}}}}),
  }};
  constexpr ElementVertexSteps<points_per_quadratic_tetrahedron> quintet_B = {{
      to_quadratic_tetrahedron(
          {{{{0, 0, 0}}, {{1, 0, 0}}, {{1, 1, 0}}, {{1, 0, 1}}}}),
      to_quadratic_tetrahedron(
          {{{{1, 1, 0}}, {{0, 1, 1}}, {{1, 0, 1}}, {{1, 1, 1}}}}),
      to_quadratic_tetrahedron(
          {{{{0, 0, 0}}, {{0, 1, 1}}, {{1, 1, 0}}, {{0, 1, 0}}}}),
      to_quadratic_tetrahedron(
          {{{{0, 0, 0}}, {{1, 0, 1}}, {{0, 1, 1}}, {{0, 0, 1}}}}),
      to_quadratic_tetrahedron(
          {{{{1, 0, 1}}, {{0, 1, 1}}, {{1, 1, 0}}, {{0, 0, 0}}}}),
  }};

  const auto add_quintet =
      [&](const ElementVertexSteps<points_per_quadratic_tetrahedron> &quintet,
          const Index &offset) {
        for (auto &&element : quintet) {
          connectivity.emplace_back();
          connectivity.back().reserve(points_per_quadratic_tetrahedron);
          for (auto &&half_steps : element) {
            connectivity.back().push_back(
                half_steps_to_id(add(half_steps, offset), granularity));
          }
        }
      };

  auto offset = Index();
  for (offset[2] = 0; offset[2] < granularity[2]; ++offset[2])
    for (offset[1] = 0; offset[1] < granularity[1]; ++offset[1])
      for (offset[0] = 0; offset[0] < granularity[0]; ++offset[0]) {
        const auto &quintet =
            (offset[0] + offset[1] + offset[2]) % 2 ? quintet_B : quintet_A;
        add_quintet(quintet, scale(offset));
      }

  return connectivity;
}

} // namespace

Mesh<Point>
generate_quadratic_tetrahedron_mesh(const Point &size,
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