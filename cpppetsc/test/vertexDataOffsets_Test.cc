// © 2021 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/cpppetsc/Mesh.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"
#include "ae108/cpppetsc/vertexDataOffsets.h"
#include <array>
#include <gmock/gmock.h>

using testing::AnyOf;
using testing::ElementsAre;
using testing::Test;

namespace ae108 {
namespace cpppetsc {
namespace {

template <class Policy> struct vertexDataOffsets_Test : Test {
  using mesh_type = Mesh<Policy>;
  using size_type = typename mesh_type::size_type;

  static constexpr size_type numberOfVertices = size_type{2};
  static constexpr size_type numberOfElements = size_type{2};
  static constexpr size_type dimension = size_type{1};
  static constexpr size_type dofPerVertex = size_type{3};
  static constexpr size_type verticesPerElement = size_type{1};

  using Connectivity =
      std::array<std::array<size_type, verticesPerElement>, numberOfElements>;
  Connectivity connectivity = {{
      {{0}},
      {{1}},
  }}; // two elements with a single vertex each
  mesh_type mesh = mesh_type::fromConnectivity(dimension, connectivity,
                                               numberOfVertices, dofPerVertex);
};

using Policies = testing::Types<SequentialComputePolicy, ParallelComputePolicy>;

TYPED_TEST_CASE(vertexDataOffsets_Test, Policies);

TYPED_TEST(vertexDataOffsets_Test,
           dofPerVertex_returns_offset_of_dofs_per_vertex) {
  const auto offsets = vertexDataOffsets(this->mesh);

  const auto dofs = this->dofPerVertex;
  EXPECT_THAT(offsets, AnyOf(ElementsAre(0, dofs), ElementsAre(dofs, 0)));
}

} // namespace
} // namespace cpppetsc
} // namespace ae108