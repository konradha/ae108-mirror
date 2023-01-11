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
#include "ae108/meshing/cppgmsh/generate_mesh.h"
#include "ae108/meshing/cppgmsh/set_granularity.h"
#include <gmock/gmock.h>
#include <gmsh.h>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct generate_mesh_Test : Test {
  void add_box() {
    gmsh::model::occ::addBox(0, 0, 0, 1, 1, 1, 1);
    gmsh::model::occ::synchronize();
  }
  void add_rectangle() {
    gmsh::model::occ::addRectangle(0, 0, 0, 1, 1, 1);
    gmsh::model::occ::synchronize();
  }

  std::size_t number_of_nodes(const int dimension) {
    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoordinates, parametricNodeCoordinates;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoordinates,
                                parametricNodeCoordinates, dimension, -1, true,
                                false);
    return nodeTags.size();
  }

  std::size_t number_of_elements(const int dimension, const int elementType) {
    std::vector<int> elementTypes{elementType};

    std::vector<std::vector<std::size_t>> elementTypeTags, elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTypeTags,
                                   elementNodeTags, dimension);

    return [&elementTypeTags]() {
      auto numberOfElements = std::size_t{0};
      for (const auto &elementTags : elementTypeTags)
        numberOfElements += elementTags.size();
      return numberOfElements;
    }();
  }
};

TEST_F(generate_mesh_Test, generate_Tet4_mesh) {

  const auto gmshContext = Context(0, 0);
  generate_mesh_Test::add_box();
  set_granularity(1.);

  const int dimension = 3;
  const int order = 1;
  const int algorithm = 6;
  generate_mesh(dimension, order, algorithm);

  const auto element_type = int(4); // 4-node tetrahedron.
  ASSERT_THAT(generate_mesh_Test::number_of_elements(dimension, element_type),
              Eq(24));
  ASSERT_THAT(generate_mesh_Test::number_of_nodes(dimension), Eq(14));
}

TEST_F(generate_mesh_Test, generate_Tet10_mesh) {

  const auto gmshContext = Context(0, 0);
  generate_mesh_Test::add_box();
  set_granularity(1.);

  const int dimension = 3;
  const int order = 2;
  const int algorithm = 6;
  generate_mesh(dimension, order, algorithm);

  const auto element_type = int(11); // 10-node second order tetrahedron
  ASSERT_THAT(generate_mesh_Test::number_of_elements(dimension, element_type),
              Eq(24));
  ASSERT_THAT(generate_mesh_Test::number_of_nodes(dimension), Eq(63));
}

TEST_F(generate_mesh_Test, generate_Tri3_mesh) {
  const auto gmshContext = Context(0, 0);

  generate_mesh_Test::add_rectangle();
  set_granularity(1.);

  const int dimension = 2;
  const int order = 1;
  const int algorithm = 6;
  generate_mesh(dimension, order, algorithm);

  const auto element_type = int(2); // 3-node triangle.
  ASSERT_THAT(generate_mesh_Test::number_of_elements(dimension, element_type),
              Eq(4));
  ASSERT_THAT(generate_mesh_Test::number_of_nodes(dimension), Eq(5));
}

TEST_F(generate_mesh_Test, generate_Tri6_mesh) {

  const auto gmshContext = Context(0, 0);
  generate_mesh_Test::add_rectangle();
  set_granularity(1.);

  const int dimension = 2;
  const int order = 2;
  const int algorithm = 6;
  generate_mesh(dimension, order, algorithm);

  const auto element_type = int(9); // 6-node second order triangle
  ASSERT_THAT(generate_mesh_Test::number_of_elements(dimension, element_type),
              Eq(4));
  ASSERT_THAT(generate_mesh_Test::number_of_nodes(dimension), Eq(13));
}

TEST_F(generate_mesh_Test, generate_Quad4_mesh) {

  const auto gmshContext = Context(0, 0);
  generate_mesh_Test::add_rectangle();
  set_granularity(1.);
  gmsh::option::setNumber("Mesh.RecombineAll", 1);
  gmsh::option::setNumber("Mesh.RecombinationAlgorithm", 2);

  const int dimension = 2;
  const int order = 1;
  const int algorithm = 6;
  generate_mesh(dimension, order, algorithm);

  const auto element_type = int(3); // 4-node quadrangle.
  ASSERT_THAT(generate_mesh_Test::number_of_elements(dimension, element_type),
              Eq(4));
  ASSERT_THAT(generate_mesh_Test::number_of_nodes(dimension), Eq(9));
}

TEST_F(generate_mesh_Test, generate_Quad9_mesh) {

  const auto gmshContext = Context(0, 0);
  generate_mesh_Test::add_rectangle();
  set_granularity(1.);
  gmsh::option::setNumber("Mesh.RecombineAll", 1);
  gmsh::option::setNumber("Mesh.RecombinationAlgorithm", 2);

  const int dimension = 2;
  const int order = 2;
  const int algorithm = 6;
  generate_mesh(dimension, order, algorithm);

  const auto element_type = int(10); // 9-node quadrangle.
  ASSERT_THAT(generate_mesh_Test::number_of_elements(dimension, element_type),
              Eq(4));
  ASSERT_THAT(generate_mesh_Test::number_of_nodes(dimension), Eq(25));
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108