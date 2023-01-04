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

#include "ae108/cpppetsc/Matrix.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cppptest/Matchers.h"
#include "ae108/solve/AffineTransform.h"
#include <gmock/gmock.h>

using ae108::cppptest::ScalarEqIfLocal;
using testing::SizeIs;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace solve {
namespace {

template <class Policy> struct AffineTransform_Test : Test {
  using transform_type = AffineTransform<Policy>;
  using vector_type = typename transform_type::vector_type;
  using matrix_type = typename transform_type::matrix_type;

  transform_type transform{
      matrix_type::fromList({{1.}, {2.}}),
      cpppetsc::tag<cpppetsc::DistributedTag>(vector_type::fromList({0., 3.}))};
};

using Policies =
    Types<cpppetsc::SequentialComputePolicy, cpppetsc::ParallelComputePolicy>;

TYPED_TEST_CASE(AffineTransform_Test, Policies);

TYPED_TEST(AffineTransform_Test, applies_transform_and_returns_vector) {
  using vector_type = typename TestFixture::vector_type;

  const auto x =
      cpppetsc::tag<cpppetsc::DistributedTag>(vector_type::fromList({-1.}));

  const auto result = apply(this->transform, x);

  ASSERT_THAT(result.unwrap(), SizeIs(2.));
  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(0, -1.));
  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(1, 1.));
}

TYPED_TEST(AffineTransform_Test,
           applies_transform_and_writes_to_output_parameter) {
  using vector_type = typename TestFixture::vector_type;

  const auto x =
      cpppetsc::tag<cpppetsc::DistributedTag>(vector_type::fromList({-1.}));

  auto result = cpppetsc::tag<cpppetsc::DistributedTag>(
      vector_type::fromLayoutOf(this->transform.vector));
  apply(this->transform, x, &result);

  ASSERT_THAT(result.unwrap(), SizeIs(2.));
  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(0, -1.));
  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(1, 1.));
}
} // namespace
} // namespace solve
} // namespace ae108