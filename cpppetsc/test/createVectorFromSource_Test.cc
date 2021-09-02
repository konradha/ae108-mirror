// © 2020 ETH Zurich, Mechanics and Materials Lab
// © 2020 California Institute of Technology
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
#include "ae108/cpppetsc/createVectorFromSource.h"
#include "ae108/cppptest/Matchers.h"
#include <gmock/gmock.h>

using ae108::cppptest::ScalarEqIfLocal;
using testing::DoubleEq;
using testing::SizeIs;
using testing::Test;

namespace ae108 {
namespace cpppetsc {
namespace {

template <class Policy> struct createVectorFromSource_Test : Test {
  using mesh_type = Mesh<Policy>;
  using size_type = typename mesh_type::size_type;
  using value_type = typename mesh_type::value_type;

  static constexpr size_type dimension = size_type{1};

  static constexpr size_type numberOfVertices = size_type{1};
  static constexpr size_type numberOfElements = size_type{1};
  static constexpr size_type verticesPerElement = size_type{1};

  static constexpr size_type dofPerVertex = size_type{0};
  static constexpr size_type dofPerElement = size_type{0};

  using Connectivity =
      std::array<std::array<size_type, verticesPerElement>, numberOfElements>;
  const Connectivity connectivity = {{{{0}}}};

  mesh_type mesh = mesh_type::fromConnectivity(
      dimension, connectivity, numberOfVertices, dofPerVertex, dofPerElement);
};

using Policies = testing::Types<SequentialComputePolicy, ParallelComputePolicy>;
TYPED_TEST_CASE(createVectorFromSource_Test, Policies);

TYPED_TEST(createVectorFromSource_Test, creates_correct_vector_with_two_dofs) {
  using size_type = typename TestFixture::size_type;
  using value_type = typename TestFixture::value_type;

  constexpr auto dofs = size_type{2};
  constexpr value_type increment = .1;

  using DataSource = std::function<void(size_type, value_type *)>;
  const auto result = createVectorFromSource(
      this->mesh, dofs,
      DataSource([&](const size_type index, value_type *const data) {
        data[0] = static_cast<value_type>(index);
        data[1] = static_cast<value_type>(index) + increment;
      }));

  ASSERT_THAT(result.unwrap(), SizeIs(dofs * this->numberOfVertices));
  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(0, 0.));
  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(1, 0. + increment));
}

} // namespace
} // namespace cpppetsc
} // namespace ae108
