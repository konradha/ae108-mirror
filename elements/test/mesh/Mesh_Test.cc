// © 2020 ETH Zurich, Mechanics and Materials Lab
//
// This file is part of ae108.
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

#include "ae108/elements/mesh/Mesh.h"
#include "ae108/elements/tensor/Tensor.h"
#include <gmock/gmock.h>
#include <stdexcept>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace elements {
namespace mesh {
namespace {

struct Mesh_Test : Test {
  using mesh_type = Mesh<tensor::Tensor<double, 2>>;
  mesh_type::Connectivity connectivity = {{
      {{0, 1, 3}},
      {{1, 2, 3}},
  }};
  mesh_type::Positions positions = {{
      {{0.0, 0.0}},
      {{1.0, 0.0}},
      {{1.0, 1.0}},
      {{0.0, 1.0}},
  }};

  mesh_type mesh{connectivity, positions};
};

TEST_F(Mesh_Test, returns_correct_connectivity) {
  EXPECT_THAT(mesh.connectivity(), Eq(connectivity));
}

TEST_F(Mesh_Test, returns_correct_position_at_index_0) {
  EXPECT_THAT(mesh.position_of_vertex(0), Eq(positions.at(0)));
}

TEST_F(Mesh_Test, returns_correct_position_at_index_1) {
  EXPECT_THAT(mesh.position_of_vertex(1), Eq(positions.at(1)));
}

TEST_F(Mesh_Test, throws_out_of_range_if_index_too_high) {
  EXPECT_THROW(mesh.position_of_vertex(123), std::out_of_range);
}

} // namespace
} // namespace mesh
} // namespace elements
} // namespace ae108