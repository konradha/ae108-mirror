// © 2021 ETH Zurich, Mechanics and Materials Lab
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