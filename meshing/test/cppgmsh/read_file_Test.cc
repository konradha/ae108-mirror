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
#include "ae108/meshing/cppgmsh/extract_mesh.h"
#include "ae108/meshing/cppgmsh/generate_mesh.h"
#include "ae108/meshing/cppgmsh/read_file.h"
#include "ae108/meshing/cppgmsh/set_granularity.h"
#include "ae108/meshing/cppgmsh/write_file.h"
#include <gmock/gmock.h>
#include <gmsh.h>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct read_file_Test : Test {};

TEST_F(read_file_Test, reads_step_file) {

  {
    const auto gmshContext = Context(0, 0);
    construct_box({0, 0, 0}, {1, 2, 3});
    gmsh::model::occ::synchronize();
    write_file("test.stp");
  }

  const auto gmshContext = Context(0, 0);

  read_file("test.stp");
  gmsh::vectorpair dimTags;
  gmsh::model::getEntities(dimTags, 3);

  ASSERT_THAT(dimTags.size(), Eq(1));
  ASSERT_THAT(dimTags[0].first, Eq(3));
  ASSERT_THAT(dimTags[0].second, Eq(1));
}

TEST_F(read_file_Test, reads_stl_file) {

  {
    const auto gmshContext = Context(0, 0);
    construct_box({0, 0, 0}, {1, 1, 1});
    gmsh::model::occ::synchronize();
    set_granularity(1.);
    generate_mesh(2, 1, 6);
    write_file("test.stl");
  }

  const auto gmshContext = Context(0, 0);

  read_file("test.stl");
  const auto [positions, connectivity, mapA, mapB] = extract_mesh<2>();

  ASSERT_THAT(positions.size(), Eq(14));
  ASSERT_THAT(connectivity.size(), Eq(24));
  ASSERT_THAT(connectivity[0].size(), Eq(3));
}

TEST_F(read_file_Test, reads_vtk_file) {

  {
    const auto gmshContext = Context(0, 0);
    construct_box({0, 0, 0}, {1, 1, 1});
    gmsh::model::occ::synchronize();
    set_granularity(1.);
    generate_mesh(3, 1, 6);
    write_file("test.vtk");
  }

  const auto gmshContext = Context(0, 0);

  read_file("test.vtk");
  gmsh::model::mesh::createTopology();
  gmsh::model::mesh::classifySurfaces(0.1, true, false, 0.1);
  gmsh::model::mesh::createGeometry();

  const auto [positions, connectivity, mapA, mapB] = extract_mesh<3>();

  ASSERT_THAT(positions.size(), Eq(14));
  ASSERT_THAT(connectivity.size(), Eq(24));
  ASSERT_THAT(connectivity[0].size(), Eq(4));
}

TEST_F(read_file_Test, reads_msh_file) {

  {
    const auto gmshContext = Context(0, 0);
    construct_box({0, 0, 0}, {1, 1, 1});
    gmsh::model::occ::synchronize();
    set_granularity(1.);
    generate_mesh(3, 1, 6);
    write_file("test.msh");
  }

  const auto gmshContext = Context(0, 0);

  read_file("test.msh");
  const auto [positions, connectivity, mapA, mapB] = extract_mesh<3>();

  ASSERT_THAT(positions.size(), Eq(14));
  ASSERT_THAT(connectivity.size(), Eq(24));
  ASSERT_THAT(connectivity[0].size(), Eq(4));
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108