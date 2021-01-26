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

#include "ae108/cpppetsc/Matrix.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cppptest/Matchers.h"
#include "ae108/solve/AffineTransform.h"
#include <gmock/gmock.h>

using ae108::cppptest::AlmostEqIfLocal;
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
  using transform_type = typename TestFixture::transform_type;
  using vector_type = typename TestFixture::vector_type;

  const auto x =
      cpppetsc::tag<cpppetsc::DistributedTag>(vector_type::fromList({-1.}));

  const auto result = apply(this->transform, x);

  ASSERT_THAT(result.unwrap(), SizeIs(2.));
  EXPECT_THAT(result.unwrap(), AlmostEqIfLocal(0, -1.));
  EXPECT_THAT(result.unwrap(), AlmostEqIfLocal(1, 1.));
}

TYPED_TEST(AffineTransform_Test,
           applies_transform_and_writes_to_output_parameter) {
  using transform_type = typename TestFixture::transform_type;
  using vector_type = typename TestFixture::vector_type;

  const auto x =
      cpppetsc::tag<cpppetsc::DistributedTag>(vector_type::fromList({-1.}));

  auto result = cpppetsc::tag<cpppetsc::DistributedTag>(
      vector_type::fromLayoutOf(this->transform.vector));
  apply(this->transform, x, &result);

  ASSERT_THAT(result.unwrap(), SizeIs(2.));
  EXPECT_THAT(result.unwrap(), AlmostEqIfLocal(0, -1.));
  EXPECT_THAT(result.unwrap(), AlmostEqIfLocal(1, 1.));
}
} // namespace
} // namespace solve
} // namespace ae108