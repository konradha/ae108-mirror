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

#include "ae108/meshing/BoundaryRepresentation.h"
#include <gmock/gmock.h>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
struct BoundaryRepresentation3D_Test : Test {

  using brep_type = BoundaryRepresentation<int, double, 3>;
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
      {{1, 3}},
      {{2, 3}},
      {{2, 1}},
  }};
  faces_type faces = {{
      {{0, 1, 5}},
      {{0, 1, 3}},
      {{0, 2, 3}},
      {{1, 2, 3}},
  }};

  brep_type brep{vertices, edges, faces};
};

struct BoundaryRepresentation2D_Test : Test {

  using brep_type = BoundaryRepresentation<int, double, 2>;
  using vertices_type = typename brep_type::Vertices;
  using edges_type = typename brep_type::Edges;
  using faces_type = typename brep_type::Faces;

  vertices_type vertices = {{
      {{0, 0}},
      {{1, 0}},
      {{0, 1}},
  }};
  edges_type edges = {{
      {{0, 1}},
      {{0, 2}},
      {{2, 1}},
  }};
  faces_type faces = {{
      {{0, 1, 2}},
  }};

  brep_type brep{vertices, edges, faces};
};

TEST_F(BoundaryRepresentation3D_Test, returns_correct_dimension) {
  EXPECT_THAT(brep.dimension(), Eq(3));
};

TEST_F(BoundaryRepresentation3D_Test, returns_correct_vertices_of_edge) {
  EXPECT_THAT(brep.vertices_of_edge(3)[0], Eq(vertices.at(1)));
  EXPECT_THAT(brep.vertices_of_edge(3)[1], Eq(vertices.at(3)));
};

TEST_F(BoundaryRepresentation3D_Test, returns_correct_vertices_of_face) {
  EXPECT_THAT(brep.vertices_of_face(0)[0], Eq(vertices.at(0)));
  EXPECT_THAT(brep.vertices_of_face(0)[1], Eq(vertices.at(1)));
  EXPECT_THAT(brep.vertices_of_face(0)[2], Eq(vertices.at(2)));
};

TEST_F(BoundaryRepresentation2D_Test, returns_correct_dimension) {
  EXPECT_THAT(brep.dimension(), Eq(2));
};

TEST_F(BoundaryRepresentation2D_Test, returns_correct_vertices_of_edge) {
  EXPECT_THAT(brep.vertices_of_edge(2)[0], Eq(vertices.at(2)));
  EXPECT_THAT(brep.vertices_of_edge(2)[1], Eq(vertices.at(1)));
};

TEST_F(BoundaryRepresentation2D_Test, returns_correct_vertices_of_face) {
  EXPECT_THAT(brep.vertices_of_face(0)[0], Eq(vertices.at(0)));
  EXPECT_THAT(brep.vertices_of_face(0)[1], Eq(vertices.at(1)));
  EXPECT_THAT(brep.vertices_of_face(0)[2], Eq(vertices.at(2)));
};

} // namespace meshing
} // namespace ae108