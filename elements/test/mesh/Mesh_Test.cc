// © 2020 ETH Zurich, Mechanics and Materials Lab
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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