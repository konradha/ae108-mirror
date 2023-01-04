// © 2020 ETH Zurich, Mechanics and Materials Lab
// © 2020 California Institute of Technology
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
