// © 2021 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/meshing/construct_periodic_point_cloud.h"
#include "ae108/meshing/construct_voronoi_cell.h"
#include <gmock/gmock.h>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
struct construct_voronoi_cell_Test : Test {};

TEST_F(construct_voronoi_cell_Test, returns_correct_cube) {

  const auto point_cloud = construct_periodic_point_cloud<3>(
      {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}, {0, 0, 0});

  std::vector<std::pair<std::size_t, std::size_t>> periodic_faces;
  const auto voronoi_cell =
      construct_voronoi_cell(point_cloud, &periodic_faces);

  ASSERT_THAT(voronoi_cell.vertices.size(), Eq(8));
  ASSERT_THAT(voronoi_cell.edges.size(), Eq(12));
  ASSERT_THAT(voronoi_cell.faces.size(), Eq(6));
  ASSERT_THAT(periodic_faces.size(), Eq(3));
};

TEST_F(construct_voronoi_cell_Test, returns_correct_hexagonal_prism) {

  const auto point_cloud = construct_periodic_point_cloud<3>(
      {{{sqrt(3.) / 2, 0.5, 0}, {sqrt(3.) / 2, -0.5, 0}, {0, 0, 1}}},
      {0, 0, 0});

  std::vector<std::pair<std::size_t, std::size_t>> periodic_faces;
  const auto voronoi_cell =
      construct_voronoi_cell(point_cloud, &periodic_faces);

  ASSERT_THAT(voronoi_cell.vertices.size(), Eq(12));
  ASSERT_THAT(voronoi_cell.edges.size(), Eq(18));
  ASSERT_THAT(voronoi_cell.faces.size(), Eq(8));
  ASSERT_THAT(periodic_faces.size(), Eq(4));
};

TEST_F(construct_voronoi_cell_Test, returns_correct_bitruncated_octahedron) {

  const auto point_cloud = construct_periodic_point_cloud<3>(
      {{{0.5, 0.5, 0.5}, {0.5, -0.5, 0.5}, {0.5, 0.5, -0.5}}}, {0, 0, 0});

  std::vector<std::pair<std::size_t, std::size_t>> periodic_faces;
  const auto voronoi_cell =
      construct_voronoi_cell(point_cloud, &periodic_faces);

  ASSERT_THAT(voronoi_cell.vertices.size(), Eq(24));
  ASSERT_THAT(voronoi_cell.edges.size(), Eq(36));
  ASSERT_THAT(voronoi_cell.faces.size(), Eq(14));
  ASSERT_THAT(periodic_faces.size(), Eq(7));
};

} // namespace meshing
} // namespace ae108