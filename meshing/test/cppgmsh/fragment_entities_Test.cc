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
#include "ae108/meshing/cppgmsh/fragment_entities.h"
#include "ae108/meshing/cppgmsh/generate_mesh.h"
#include "ae108/meshing/cppgmsh/get_boundary_of.h"
#include "ae108/meshing/cppgmsh/get_elements_in.h"
#include "ae108/meshing/cppgmsh/get_nodes_of.h"
#include "ae108/meshing/cppgmsh/intersect_entities.h"
#include "ae108/meshing/cppgmsh/set_granularity.h"
#include "ae108/meshing/cppgmsh/synchronize.h"
#include <gmock/gmock.h>
#include <gmsh.h>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct fragment_Test : Test {};

TEST_F(fragment_Test, returns_the_correct_entities) {
  const auto gmshContext = Context(0, 0);

  const auto box_1 = construct_box({0., 0., 0.}, {1., 1., 1.});
  const auto box_2 = construct_box({0., 0., 1.}, {1., 1., 1.});
  const auto fragmented_entities = fragment_entities({box_1, box_2});

  ASSERT_THAT(fragmented_entities.size(), Eq(2));
  ASSERT_THAT(fragmented_entities[0], Eq(box_1));
  ASSERT_THAT(fragmented_entities[1], Eq(box_2));
}

TEST_F(fragment_Test, no_elements_crossing_entity_interface) {
  const auto gmshContext = Context(0, 0);

  const auto box_1 = construct_box({0., 0., 0.}, {1., 1., 6. / 7.});
  const auto box_2 = construct_box({0., 0., 6. / 7.}, {1., 1., 1. / 7.});
  const auto fragmented_entities = fragment_entities({box_1, box_2});
  synchronize();

  set_granularity(0.5);
  generate_mesh(3, 1, 6);

  const auto elements = get_elements_in({3, -1});
  const auto elements_1 = get_elements_in(box_1);
  const auto elements_2 = get_elements_in(box_2);
  ASSERT_THAT(elements_1.size() + elements_2.size(), Eq(elements.size()));

  const auto nodes_1 = get_nodes_of<3>(box_1);
  for (const auto &element_1 : elements_1) {

    int elementType;
    std::vector<std::size_t> nodes;
    gmsh::model::mesh::getElement(element_1.second, elementType, nodes);

    bool is_in_nodes_1 = false;
    for (const auto &node_1 : nodes_1)
      is_in_nodes_1 += std::count(nodes.begin(), nodes.end(), node_1.id);

    ASSERT_THAT(is_in_nodes_1, Eq(true));
  }

  auto nodes_2 = get_nodes_of<3>(box_2);
  for (const auto &element_2 : elements_2) {

    int elementType;
    std::vector<std::size_t> nodes;
    gmsh::model::mesh::getElement(element_2.second, elementType, nodes);

    bool is_in_nodes_2 = false;
    for (const auto &node_2 : nodes_2)
      is_in_nodes_2 += std::count(nodes.begin(), nodes.end(), node_2.id);

    ASSERT_THAT(is_in_nodes_2, Eq(true));
  }
}

TEST_F(fragment_Test, conformal_mesh_at_interface) {
  const auto gmshContext = Context(0, 0);

  const auto box_1 = construct_box({0., 0., 0.}, {1., 1., 6. / 7.});
  const auto box_2 =
      construct_box({1. / 3., 1. / 3., 6. / 7.}, {1. / 3., 1. / 3., 1. / 7.});
  synchronize();

  const auto boundary_1 = get_boundary_of({box_1});
  const auto boundary_2 = get_boundary_of({box_2});

  const auto common_boundary = intersect_entities(boundary_1, boundary_2);
  const auto fragmented_entities = fragment_entities({box_1, box_2});
  synchronize();

  set_granularity(0.5);
  generate_mesh(3, 1, 6);

  const auto nodes_1 = get_nodes_of<3>(box_1);
  const auto nodes_2 = get_nodes_of<3>(box_2);
  const auto nodes_interface = get_nodes_of<3>(common_boundary[0]);

  bool is_in_nodes_1 = false;
  for (const auto &node_1 : nodes_1)
    is_in_nodes_1 +=
        std::count(nodes_interface.begin(), nodes_interface.end(), node_1);
  ASSERT_THAT(is_in_nodes_1, Eq(true));

  bool is_in_nodes_2 = false;
  for (const auto &node_2 : nodes_2)
    is_in_nodes_2 +=
        std::count(nodes_interface.begin(), nodes_interface.end(), node_2);
  ASSERT_THAT(is_in_nodes_2, Eq(true));
}

TEST_F(fragment_Test, handles_overlapping_entities) {
  const auto gmshContext = Context(0, 0);

  const auto box_1 = construct_box({0., 0., 0.}, {1., 1., 1.});
  const auto box_2 = construct_box({0.5, 0.5, 0.5}, {1., 1., 1.});
  const auto fragmented_entities = fragment_entities({box_1, box_2});

  ASSERT_THAT(fragmented_entities.size(), Eq(3));
}

TEST_F(fragment_Test, handles_unconnected_entities) {
  const auto gmshContext = Context(0, 0);

  const auto box_1 = construct_box({0., 0., 0.}, {1., 1., 1.});
  const auto box_2 = construct_box({2., 2., 2.}, {1., 1., 1.});
  auto fragmented_entities = fragment_entities({box_1, box_2});

  ASSERT_THAT(fragmented_entities.size(), Eq(2));
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108