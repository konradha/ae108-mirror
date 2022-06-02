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

#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cpppetsc/nestVectors.h"
#include "ae108/cppptest/Matchers.h"
#include <cmath>
#include <gmock/gmock.h>

using testing::DoubleEq;
using testing::Eq;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace cpppetsc {

template <class Policy> struct nestVectors_Test : Test {
  using vector_type = Vector<Policy>;
};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;

TYPED_TEST_CASE(nestVectors_Test, Policies);

TYPED_TEST(nestVectors_Test, result_has_sum_of_rows) {
  using vector_type = typename TestFixture::vector_type;

  auto x = distributed<vector_type>(2);
  auto y = distributed<vector_type>(3);

  const auto nested = nestVectors<TypeParam>({&x, &y});

  EXPECT_THAT(nested.unwrap().size(), Eq(2 + 3));
}

TYPED_TEST(nestVectors_Test, result_has_correct_magnitude) {
  using vector_type = typename TestFixture::vector_type;

  auto x = distributed<vector_type>(vector_type::fromList({1., 2.}));
  auto y = distributed<vector_type>(vector_type::fromList({3., 4., 5.}));

  const auto nested = nestVectors<TypeParam>({&x, &y});

  EXPECT_THAT(std::pow(nested.unwrap().norm(), 2.),
              DoubleEq(std::pow(x.unwrap().norm(), 2.) +
                       std::pow(y.unwrap().norm(), 2.)));
}

} // namespace cpppetsc
} // namespace ae108