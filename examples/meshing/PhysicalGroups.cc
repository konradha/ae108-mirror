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
#include "ae108/meshing/cppgmsh/get_entities_in.h"
#include "ae108/meshing/cppgmsh/get_nodes_of.h"
#include "ae108/meshing/cppgmsh/get_physical_groups.h"
#include "ae108/meshing/cppgmsh/set_granularity.h"
#include "ae108/meshing/cppgmsh/set_physical_group_of.h"
#include "ae108/meshing/cppgmsh/synchronize.h"

#include <array>
#include <iostream>

using namespace ae108::meshing;

/*
In this example we will construct and mesh a cube with unit side lengths, which
consists of three layers (z < 0.3 -> layer 1, 0.3 < z < 0.7 -> layer 2, z > 0.7
-> layer 3). Layers 1 and 3 are of one material, layer 2 of another. In order to
account for this, physical groups are defined.

     ________________    _
    |                |  |
    |   Material 1   |  | 0.3
    |________________|  |_
    |                |  |
    |   Material 2   |  | 0.4
    |                |  |
    |________________|  |_
    |                |  |
    |   Material 1   |  | 0.3
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

  // Construct three boxes:
  const auto boxes = std::vector<std::pair<int, int>>{
      cppgmsh::construct_box({0., 0., 0.}, {1., 1., 0.3}),
      cppgmsh::construct_box({0., 0., .3}, {1., 1., 0.4}),
      cppgmsh::construct_box({0., 0., .7}, {1., 1., 0.3})};

  std::cout << "Box 1: dim = " << boxes[0].first
            << ", tag = " << boxes[0].second << std::endl
            << "Box 2: dim = " << boxes[1].first
            << ", tag = " << boxes[1].second << std::endl
            << "Box 3: dim = " << boxes[2].first
            << ", tag = " << boxes[2].second << std::endl;

  // Fuse boxes, while keeping track of the original entities:
  std::vector<std::pair<int, int>> fused_entitites =
      cppgmsh::fragment_entities(boxes);
  cppgmsh::synchronize();

  // Define physical groups:
  cppgmsh::set_physical_group_of({boxes[0], boxes[2]});
  cppgmsh::set_physical_group_of({boxes[1]});

  // Show that the physical groups have been defined as desired:
  auto groups = cppgmsh::get_physical_groups();
  std::cout << "Groups (dim, tag): ";
  for (const auto &group : groups)
    std::cout << "(" << group.first << "," << group.second << ") ";
  std::cout << std::endl;

  // Get the tags of the entities attributed two both groups:
  for (const auto &group : groups) {
    std::cout << "Group " << group.second << ", dim: " << group.first
              << ", entities: ";
    const auto entities = cppgmsh::get_entities_in(group);
    for (const auto &entity : entities)
      std::cout << entity.second << " ";
    std::cout << std::endl;
  }

  // Create mesh:
  cppgmsh::set_granularity(2.);
  cppgmsh::generate_mesh(3, 1, 6);

  const auto [positions, connectivity, nodeTagToIndex, elementTagToIndex] =
      cppgmsh::extract_mesh<3>();

  // Print element tags and indides, as well as node indides for the two groups:
  for (const auto &group : groups) {
    std::cout << "Elements in group " << group.second << ":\n";
    std::cout << "\ttag\tindex\tnode indides\n";
    const auto entities = cppgmsh::get_entities_in(group);
    for (const auto &entity : entities) {
      const auto elements = cppgmsh::get_elements_in(entity);
      for (const auto &element : elements) {
        const auto index = elementTagToIndex.at(element.second);
        std::cout << "\t" << element.second << "\t" << index << "\t";
        for (const auto &node : connectivity[index])
          std::cout << node << "\t";
        std::cout << std::endl;
      }
    }
  }

  // Print node tags, indides and coordinates for the two boxes:
  for (const auto &group : groups) {
    std::cout << "Nodes in group " << group.second << ":\n";
    std::cout << "\ttag\tindex\tx\ty\tz\n";
    const auto entities = cppgmsh::get_entities_in(group);
    for (const auto &entity : entities) {
      const auto nodes = cppgmsh::get_nodes_of<3>(entity);
      for (const auto &node : nodes) {
        std::cout << "\t" << node.id << "\t" << nodeTagToIndex.at(node.id)
                  << "\t" << node.position[0] << "\t" << node.position[1]
                  << "\t" << node.position[2] << std::endl;
      }
    }
  }
}