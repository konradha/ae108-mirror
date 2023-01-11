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

#include "ae108/meshing/cppgmsh/Context.h"
#include "ae108/meshing/cppgmsh/construct_box.h"
#include "ae108/meshing/cppgmsh/construct_rectangle.h"
#include "ae108/meshing/cppgmsh/extract_mesh.h"
#include "ae108/meshing/cppgmsh/generate_mesh.h"
#include "ae108/meshing/cppgmsh/set_granularity.h"
#include "ae108/meshing/cppgmsh/synchronize.h"
#include <gmock/gmock.h>
#include <gmsh.h>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct extract_mesh_Test : Test {};

TEST_F(extract_mesh_Test, extract_Tet4_mesh) {
  const auto gmshContext = Context(0, 0);
  construct_box({0, 0, 0}, {1, 1, 1});
  synchronize();
  set_granularity(1.);

  const int dimension = 3;
  const int order = 1;
  const int algorithm = 6;
  generate_mesh(dimension, order, algorithm);

  const auto [positions, connectivity, mapA, mapB] = extract_mesh<dimension>();

  ASSERT_THAT(positions.size(), Eq(14));
  ASSERT_THAT(connectivity.size(), Eq(24));
  ASSERT_THAT(connectivity[0].size(), Eq(4));
}

TEST_F(extract_mesh_Test, extract_Tet10_mesh) {
  const auto gmshContext = Context(0, 0);
  construct_box({0, 0, 0}, {1, 1, 1});
  synchronize();
  set_granularity(1.);

  const int dimension = 3;
  const int order = 2;
  const int algorithm = 6;
  generate_mesh(dimension, order, algorithm);
  const auto [positions, connectivity, mapA, mapB] = extract_mesh<dimension>();

  ASSERT_THAT(positions.size(), Eq(63));
  ASSERT_THAT(connectivity.size(), Eq(24));
  ASSERT_THAT(connectivity[0].size(), Eq(10));
}

TEST_F(extract_mesh_Test, extract_Tri3_mesh) {
  const auto gmshContext = Context(0, 0);
  construct_rectangle({0, 0, 0}, {1, 1});
  synchronize();
  set_granularity(1.);

  const int dimension = 2;
  const int order = 1;
  const int algorithm = 6;
  generate_mesh(dimension, order, algorithm);
  const auto [positions, connectivity, mapA, mapB] = extract_mesh<dimension>();

  ASSERT_THAT(positions.size(), Eq(5));
  ASSERT_THAT(connectivity.size(), Eq(4));
  ASSERT_THAT(connectivity[0].size(), Eq(3));
}

TEST_F(extract_mesh_Test, extract_Tri6_mesh) {
  const auto gmshContext = Context(0, 0);
  construct_rectangle({0, 0, 0}, {1, 1});
  synchronize();
  set_granularity(1.);

  const int dimension = 2;
  const int order = 2;
  const int algorithm = 6;
  generate_mesh(dimension, order, algorithm);
  const auto [positions, connectivity, mapA, mapB] = extract_mesh<dimension>();

  ASSERT_THAT(positions.size(), Eq(13));
  ASSERT_THAT(connectivity.size(), Eq(4));
  ASSERT_THAT(connectivity[0].size(), Eq(6));
}

TEST_F(extract_mesh_Test, extract_Quad4_mesh) {
  const auto gmshContext = Context(0, 0);
  construct_rectangle({0, 0, 0}, {1, 1});
  synchronize();
  set_granularity(1.);
  gmsh::option::setNumber("Mesh.RecombineAll", 1);
  gmsh::option::setNumber("Mesh.RecombinationAlgorithm", 2);

  const int dimension = 2;
  const int order = 1;
  const int algorithm = 6;
  generate_mesh(dimension, order, algorithm);
  const auto [positions, connectivity, mapA, mapB] = extract_mesh<dimension>();

  ASSERT_THAT(positions.size(), Eq(9));
  ASSERT_THAT(connectivity.size(), Eq(4));
  ASSERT_THAT(connectivity[0].size(), Eq(4));
}

TEST_F(extract_mesh_Test, extract_Quad9_mesh) {
  const auto gmshContext = Context(0, 0);
  construct_rectangle({0, 0, 0}, {1, 1});
  synchronize();
  set_granularity(1.);
  gmsh::option::setNumber("Mesh.RecombineAll", 1);
  gmsh::option::setNumber("Mesh.RecombinationAlgorithm", 2);

  const int dimension = 2;
  const int order = 2;
  const int algorithm = 6;
  generate_mesh(dimension, order, algorithm);
  const auto [positions, connectivity, mapA, mapB] = extract_mesh<dimension>();

  ASSERT_THAT(positions.size(), Eq(25));
  ASSERT_THAT(connectivity.size(), Eq(4));
  ASSERT_THAT(connectivity[0].size(), Eq(9));
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108