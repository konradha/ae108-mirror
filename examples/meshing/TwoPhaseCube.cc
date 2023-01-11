// © 2022 ETH Zurich, Mechanics and Materials Lab
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
#include "ae108/meshing/cppgmsh/fragment_entities.h"
#include "ae108/meshing/cppgmsh/generate_mesh.h"
#include "ae108/meshing/cppgmsh/get_elements_in.h"
#include "ae108/meshing/cppgmsh/get_nodes_of.h"
#include "ae108/meshing/cppgmsh/set_granularity.h"
#include "ae108/meshing/cppgmsh/synchronize.h"

#include <array>
#include <iostream>

using namespace ae108::meshing;

/*
In this example we will construct and mesh a cube with unit side lengths, which
consists of two different phases (z < 0.5 -> phase 1,z > 0.5 -> phase 2).
     ________________    _
    |                |  |
    |     Phase 2    |  | 0.5
    |________________|  |_
    |                |  |
    |     Phase 1    |  | 0.5
    |________________|  |_
     ________________
    |       1.0      |
 z
 ^
 |
    --> x,y

*/

int main(int argc, char **argv) {

  // Initialize the gmsh context:
  const auto gmshContext = cppgmsh::Context(argc, argv, true, false);

  // Construct two boxes:
  const auto boxes = std::vector<std::pair<int, int>>{
      cppgmsh::construct_box({0., 0., 0.}, {1., 1., 0.5}),
      cppgmsh::construct_box({0., 0., .5}, {1., 1., 0.5})};

  std::cout << "Box 1: dim = " << boxes[0].first
            << ", tag = " << boxes[0].second << std::endl
            << "Box 2: dim = " << boxes[1].first
            << ", tag = " << boxes[1].second << std::endl;

  // Fuse boxes, while keeping track of the original entities:
  const auto fused_entitites = cppgmsh::fragment_entities(boxes);
  cppgmsh::synchronize();

  // Create mesh:
  cppgmsh::set_granularity(2.);
  cppgmsh::generate_mesh(3, 1, 6);

  const auto [positions, connectivity, nodeTagToIndex, elementTagToIndex] =
      cppgmsh::extract_mesh<3>();

  // Print element tags and indides, as well as node indides for the two boxes:
  for (const auto &box : boxes) {
    std::cout << "Elements in box " << box.second << ":\n"
              << "\ttag\tindex\tnode indices\n";

    const auto elements = cppgmsh::get_elements_in(box);
    for (const auto &element : elements) {
      const auto index = elementTagToIndex.at(element.second);
      std::cout << "\t" << element.second << "\t" << index << "\t";
      for (const auto &node : connectivity[index])
        std::cout << node << "\t";
      std::cout << std::endl;
    }
  }

  // Print node tags, indides and coordinates for the two boxes:
  for (const auto &box : boxes) {
    std::cout << "Nodes in box " << box.second << ":\n"
              << "\ttag\tindex\tx\ty\tz\n";

    const auto nodes = cppgmsh::get_nodes_of<3>(box);
    for (const auto &node : nodes)
      std::cout << "\t" << node.id << "\t" << nodeTagToIndex.at(node.id) << "\t"
                << node.position[0] << "\t" << node.position[1] << "\t"
                << node.position[2] << std::endl;
  }
}