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
#include "ae108/meshing/cppgmsh/construct_box.h"
#include "ae108/meshing/cppgmsh/construct_rectangle.h"
#include "ae108/meshing/cppgmsh/extract_mesh.h"
#include "ae108/meshing/cppgmsh/extrude_surface.h"
#include "ae108/meshing/cppgmsh/fuse_entities.h"
#include "ae108/meshing/cppgmsh/generate_mesh.h"
#include "ae108/meshing/cppgmsh/remove_entities.h"
#include "ae108/meshing/cppgmsh/rotate_entities.h"
#include "ae108/meshing/cppgmsh/set_granularity.h"
#include "ae108/meshing/cppgmsh/synchronize.h"
#include "ae108/meshing/cppgmsh/translate_entities.h"
#include "ae108/meshing/cppgmsh/write_file.h"

#include <gmsh.h>
#include <iostream>
#include <numeric>

using namespace ae108::meshing;

constexpr std::size_t dimension = 3;
using Point = std::array<double, dimension>;

int main(int argc, char **argv) {

  const auto gmshContext = cppgmsh::Context(argc, argv, true, false);

  const double gamma_start = 5.;
  const double gamma_end = 10.;
  const double length = 466.75;

  const auto xcoords = [] {
    std::vector<double> range(80);
    std::iota(range.begin(), range.end(), -40);
    std::transform(range.begin(), range.end(), range.begin(),
                   [](auto &N) { return N * 5; });
    return range;
  }();

  const auto ycoords = [&] {
    std::vector<double> range(59);
    std::iota(range.begin(), range.end(), -29);
    std::transform(range.begin(), range.end(), range.begin(), [&](auto &N) {
      return N * length * gamma_start /
             (length - fabs(N) * gamma_end + fabs(N) * gamma_start);
    });
    return range;
  }();

  const auto crosssection =
      cppgmsh::construct_rectangle({-.25, -.5, 0}, {.5, 1.});
  cppgmsh::rotate_entities({crosssection}, {0, 0, 0}, {1, 0, 0}, M_PI_2);

  std::vector<std::pair<int, int>> entities;
  for (const auto &xcoord : xcoords)
    entities.push_back(cppgmsh::extrude_surface(
        crosssection.second, {xcoord, -200, 0}, {xcoord, 200, 0}, {0, 1, 0}));
  for (const auto &ycoord : ycoords)
    entities.push_back(cppgmsh::extrude_surface(
        crosssection.second, {-200, ycoord, 0}, {200, ycoord, 0}, {0, 1, 0}));

  entities.push_back(cppgmsh::construct_box({-225, -225, -.5}, {450, 25, 1}));
  entities.push_back(cppgmsh::construct_box({-225, -225, -.5}, {25, 450, 1}));
  entities.push_back(cppgmsh::construct_box({-225, 200, -.5}, {450, 25, 1}));
  entities.push_back(cppgmsh::construct_box({200, -225, -.5}, {25, 450, 1}));

  cppgmsh::remove_entities({crosssection}, true);

  cppgmsh::fuse_entities(entities);

  cppgmsh::synchronize();
  cppgmsh::write_file("gradedLattice.step");

  cppgmsh::generate_mesh();
  cppgmsh::write_file("gradedLattice.vtk");

  const auto [positions, connectivity, _, __] =
      cppgmsh::extract_mesh<dimension>();

  std::cout << "# points: " << positions.size() << std::endl
            << "# elements: " << connectivity.size() << std::endl;
}