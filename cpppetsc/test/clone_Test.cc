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
#include "ae108/cpppetsc/TaggedVector.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cpppetsc/clone.h"
#include "ae108/cppptest/Matchers.h"
#include <gmock/gmock.h>

using ae108::cppptest::ScalarEqIfLocal;
using testing::Pair;
using testing::SizeIs;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace cpppetsc {

template <class Policy> struct clone_Test : Test {
  using vector_type = Vector<Policy>;
  using matrix_type = Matrix<Policy>;
};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;

TYPED_TEST_CASE(clone_Test, Policies);

TYPED_TEST(clone_Test, clones_vectors) {
  using vector_type = typename TestFixture::vector_type;

  const auto from = tag<DistributedTag>(vector_type::fromList({3., 7.}));

  const auto result = clone(from);

  ASSERT_THAT(result.unwrap(), SizeIs(2));
  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(0, 3.));
  EXPECT_THAT(result.unwrap(), ScalarEqIfLocal(1, 7.));
}

TYPED_TEST(clone_Test, clones_matrices) {
  using matrix_type = typename TestFixture::matrix_type;

  const auto from = matrix_type::fromList({{1., 2.}, {3., 4.}});

  const auto result = clone(from);

  ASSERT_THAT(result.size(), Pair(2, 2));
  EXPECT_THAT(result, ScalarEqIfLocal(0, 0, 1.));
  EXPECT_THAT(result, ScalarEqIfLocal(0, 1, 2.));
  EXPECT_THAT(result, ScalarEqIfLocal(1, 0, 3.));
  EXPECT_THAT(result, ScalarEqIfLocal(1, 1, 4.));
}

} // namespace cpppetsc
} // namespace ae108