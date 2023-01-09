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

#include "ae108/meshing/construct_rectilinear_grid.h"
#include "ae108/meshing/cppgmsh/Context.h"
#include "ae108/meshing/cppgmsh/construct_rectangle.h"
#include "ae108/meshing/cppgmsh/extract_mesh.h"
#include "ae108/meshing/cppgmsh/extrude_surface.h"
#include "ae108/meshing/cppgmsh/fuse_entities.h"
#include "ae108/meshing/cppgmsh/generate_mesh.h"
#include "ae108/meshing/cppgmsh/remove_entities.h"
#include "ae108/meshing/cppgmsh/set_granularity.h"
#include "ae108/meshing/cppgmsh/synchronize.h"
#include "ae108/meshing/cppgmsh/write_file.h"

#include <gmsh.h>
#include <iostream>

using namespace ae108::meshing;

constexpr std::size_t dimension = 3;
using Point = std::array<double, dimension>;

int main(int argc, char **argv) {

  const auto gmshContext = cppgmsh::Context(argc, argv, true, false);

  // construct a regular grid
  //    *--*
  //    |  |
  //    *--*

  const auto coordinates =
      std::array<std::vector<double>, dimension>{{{0, 1}, {0, 1}, {0}}};
  auto [nodes, segments] = construct_rectilinear_grid(coordinates);

  const auto crosssection =
      cppgmsh::construct_rectangle({-.05, -.05, 0}, {.1, .1});

  std::vector<std::pair<int, int>> beams;
  for (const auto &segment : segments)
    beams.push_back(cppgmsh::extrude_surface(
        crosssection.second, nodes[segment[0]], nodes[segment[1]]));

  cppgmsh::remove_entities({crosssection}, true);

  cppgmsh::fuse_entities(beams);

  cppgmsh::synchronize();
  //   cppgmsh::write_file("segmentMesh.step");

  cppgmsh::set_granularity(0.25);
  cppgmsh::generate_mesh();
  cppgmsh::write_file("segmentMesh.vtk");

  const auto [positions, connectivity, _, __] =
      cppgmsh::extract_mesh<dimension>();

  std::cout << "# points: " << positions.size() << std::endl
            << "# elements: " << connectivity.size() << std::endl;
}