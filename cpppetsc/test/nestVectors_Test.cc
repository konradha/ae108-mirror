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