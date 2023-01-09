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
#include "ae108/meshing/cppgmsh/construct_polyhedron.h"
#include "ae108/meshing/cppgmsh/get_entities_of.h"
#include "ae108/meshing/cppgmsh/synchronize.h"
#include <gmock/gmock.h>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct construct_polyhedron_Test : Test {

  using brep_type = BoundaryRepresentation<std::size_t, double, 3>;
  using vertices_type = typename brep_type::Vertices;
  using edges_type = typename brep_type::Edges;
  using faces_type = typename brep_type::Faces;

  vertices_type vertices = {{
      {{0, 0, 0}},
      {{1, 0, 0}},
      {{0, 1, 0}},
      {{0, 0, 1}},
  }};
  edges_type edges = {{
      {{0, 1}},
      {{0, 2}},
      {{0, 3}},
      {{1, 2}},
      {{1, 3}},
      {{2, 3}},
  }};
  faces_type faces = {{
      {{0, 1, 3}},
      {{0, 2, 4}},
      {{1, 2, 5}},
      {{3, 4, 5}},
  }};

  brep_type brep{vertices, edges, faces};
};

TEST_F(construct_polyhedron_Test, constructs_tetrahedron) {
  const auto gmshContext = Context(0, 0, true, true);

  const auto polyhedron = construct_polyhedron(brep);

  ASSERT_THAT(polyhedron.first, Eq(3));
  ASSERT_THAT(polyhedron.second, Eq(1));
}

TEST_F(construct_polyhedron_Test, returns_correct_number_of_volumes) {
  const auto gmshContext = Context(0, 0, true, true);

  construct_polyhedron(brep);
  synchronize();

  const auto surfaces = get_entities_of(3);

  ASSERT_THAT(surfaces.size(), Eq(1));
}

TEST_F(construct_polyhedron_Test, returns_correct_number_of_surfaces) {
  const auto gmshContext = Context(0, 0, true, true);

  construct_polyhedron(brep);
  synchronize();

  const auto surfaces = get_entities_of(2);

  ASSERT_THAT(surfaces.size(), Eq(4));
}

TEST_F(construct_polyhedron_Test, returns_correct_number_of_lines) {
  const auto gmshContext = Context(0, 0, true, true);

  construct_polyhedron(brep);
  synchronize();

  const auto surfaces = get_entities_of(1);

  ASSERT_THAT(surfaces.size(), Eq(6));
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108