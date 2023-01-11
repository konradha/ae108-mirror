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

#include "ae108/meshing/BoundingBox.h"
#include "ae108/meshing/cppgmsh/Context.h"
#include "ae108/meshing/cppgmsh/construct_box.h"
#include "ae108/meshing/cppgmsh/get_entities_in.h"
#include "ae108/meshing/cppgmsh/get_points_of.h"
#include "ae108/meshing/cppgmsh/heal_periodic_domains.h"
#include <gmock/gmock.h>
#include <gmsh.h>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct heal_periodic_domains_Test : Test {};

TEST_F(heal_periodic_domains_Test, heal_periodic_box_surfaces) {
  const auto gmshContext = Context(0, 0);
  const auto box = construct_box({0, 0, 0}, {1, 1, 1});

  gmsh::vectorpair ov;
  std::vector<gmsh::vectorpair> ovv;
  gmsh::model::occ::fragment(
      {box}, {{0, gmsh::model::occ::addPoint(0.5, 0, 0)}}, ov, ovv);
  gmsh::model::occ::synchronize();

  const auto source_bbox =
      BoundingBox<std::array<double, 3>>{{0, 0, 0}, {1, 1, 0}};
  const auto target_bbox =
      BoundingBox<std::array<double, 3>>{{0, 0, 1}, {1, 1, 1}};

  heal_periodic_domains(source_bbox, target_bbox);
  gmsh::model::occ::synchronize();

  const auto source_surface = get_entities_in(source_bbox, 2);
  const auto target_surface = get_entities_in(target_bbox, 2);
  const auto source_points = get_points_of(source_surface[0]);
  const auto target_points = get_points_of(target_surface[0]);

  gmsh::vectorpair points;
  gmsh::model::getEntities(points, 0);

  ASSERT_THAT(points.size(), Eq(10));
  ASSERT_THAT(source_points.size(), Eq(5));
  ASSERT_THAT(target_points.size(), Eq(5));
}

TEST_F(heal_periodic_domains_Test, heal_two_periodic_box_surfaces) {
  const auto gmshContext = Context(0, 0);
  const auto boxA = construct_box({0, 0, 0}, {.4, .4, 1});
  construct_box({.5, .5, 0}, {.4, .4, 1});

  gmsh::vectorpair ov;
  std::vector<gmsh::vectorpair> ovv;
  gmsh::model::occ::fragment(
      {boxA}, {{0, gmsh::model::occ::addPoint(0.25, 0, 0)}}, ov, ovv);
  gmsh::model::occ::synchronize();

  const auto source_bbox =
      BoundingBox<std::array<double, 3>>{{0, 0, 0}, {1, 1, 0}};
  const auto target_bbox =
      BoundingBox<std::array<double, 3>>{{0, 0, 1}, {1, 1, 1}};

  heal_periodic_domains(source_bbox, target_bbox);
  gmsh::model::occ::synchronize();

  const auto source_surface = get_entities_in(source_bbox, 2);
  const auto target_surface = get_entities_in(target_bbox, 2);
  const auto source_points = get_points_of(source_surface[0]).size() +
                             get_points_of(source_surface[1]).size();
  const auto target_points = get_points_of(target_surface[0]).size() +
                             get_points_of(target_surface[1]).size();

  gmsh::vectorpair points;
  gmsh::model::getEntities(points, 0);

  ASSERT_THAT(points.size(), Eq(18));
  ASSERT_THAT(source_points, Eq(9));
  ASSERT_THAT(target_points, Eq(9));
}
} // namespace cppgmsh
} // namespace meshing
} // namespace ae108