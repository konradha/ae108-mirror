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
#include "ae108/meshing/cppgmsh/generate_mesh.h"
#include "ae108/meshing/cppgmsh/get_elements_in.h"
#include "ae108/meshing/cppgmsh/get_nodes_of.h"
#include "ae108/meshing/cppgmsh/get_points_of.h"
#include "ae108/meshing/cppgmsh/set_granularity.h"
#include <gmock/gmock.h>
#include <gmsh.h>

using testing::ElementsAre;
using testing::Eq;
using testing::Ge;
using testing::Le;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct get_elements_in_Test : Test {};

TEST_F(get_elements_in_Test, returns_correct_number_of_tets_for_simple_box) {
  const auto gmshContext = Context(0, 0);

  construct_box({0., 0., 0.}, {1., 1., 1.});
  gmsh::model::occ::synchronize();

  set_granularity(1);
  generate_mesh(3, 1, 6);

  auto elements = get_elements_in({3, -1});
  ASSERT_THAT(elements.size(), Eq(24));
  ASSERT_THAT(elements[0].first, Eq(4));
}

TEST_F(get_elements_in_Test, returns_the_correct_elements) {
  const auto gmshContext = Context(0, 0);

  const auto box_1 = construct_box({0., 0., 0.}, {1., 1., 1.});
  const auto box_2 = construct_box({10., 10., 10.}, {1., 1., 1.});
  gmsh::model::occ::synchronize();

  set_granularity(1.);
  generate_mesh(3, 1, 6);

  auto elements = get_elements_in({3, -1});
  std::vector<bool> foundElement(elements.size(), false);

  auto elements_1 = get_elements_in({3, box_1.second});
  auto elements_2 = get_elements_in({3, box_2.second});

  const auto all_nodes = get_nodes_of<3>({3, -1});
  for (std::size_t i = 0; i < elements_1.size(); i++) {
    int type;
    std::vector<std::size_t> nodes;
    gmsh::model::mesh::getElement(elements_1[i].second, type, nodes);
    ASSERT_THAT(elements_1[i].first, Eq(type));
    for (std::size_t j = 0; j < nodes.size(); j++) {
      bool foundNode = false;
      for (std::size_t k = 0; k < all_nodes.size(); k++) {
        if (nodes[j] == all_nodes[k].id) {
          foundNode = true;
          ASSERT_THAT(all_nodes[k].position,
                      ElementsAre(Ge(0.), Ge(0.), Ge(0.)));
          ASSERT_THAT(all_nodes[k].position,
                      ElementsAre(Le(1.), Le(1.), Le(1.)));
          break;
        }
      }
      ASSERT_THAT(foundNode, Eq(true));
    }
    for (std::size_t j = 0; j < elements.size(); j++) {
      if (elements_1[i].second == elements[j].second) {
        ASSERT_THAT(foundElement[j], Eq(false));
        foundElement[j] = true;
      }
    }
  }

  for (std::size_t i = 0; i < elements_2.size(); i++) {
    int type;
    std::vector<std::size_t> nodes;
    gmsh::model::mesh::getElement(elements_2[i].second, type, nodes);
    for (std::size_t j = 0; j < nodes.size(); j++) {
      bool foundNode = false;
      for (std::size_t k = 0; k < all_nodes.size(); k++) {
        if (nodes[j] == all_nodes[k].id) {
          foundNode = true;
          ASSERT_THAT(all_nodes[k].position,
                      ElementsAre(Ge(10.), Ge(10.), Ge(10.)));
          ASSERT_THAT(all_nodes[k].position,
                      ElementsAre(Le(11.), Le(11.), Le(11.)));
          break;
        }
      }
      ASSERT_THAT(foundNode, Eq(true));
    }
    for (std::size_t j = 0; j < elements.size(); j++) {
      if (elements_2[i].second == elements[j].second) {
        ASSERT_THAT(foundElement[j], Eq(false));
        foundElement[j] = true;
      }
    }
  }
  for (std::size_t i = 0; i < elements.size(); i++) {
    ASSERT_THAT(foundElement[i], Eq(true));
  }
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108