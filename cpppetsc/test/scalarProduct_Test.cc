// © 2022 ETH Zurich, Mechanics and Materials Lab
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
#include "ae108/cpppetsc/scalarProduct.h"
#include "ae108/cppptest/Matchers.h"
#include <gmock/gmock.h>

using ae108::cppptest::ScalarEq;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace cpppetsc {

template <class Policy> struct scalarProduct_Test : Test {
  using vector_type = Vector<Policy>;
};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;

TYPED_TEST_CASE(scalarProduct_Test, Policies);

TYPED_TEST(scalarProduct_Test, computes_scalar_product_of_vectors) {
  using vector_type = typename TestFixture::vector_type;

  const auto x = tag<DistributedTag>(vector_type::fromList({3., 7.}));
  const auto y = tag<DistributedTag>(vector_type::fromList({2., -4.}));

  EXPECT_THAT(scalarProduct(x, y), ScalarEq(6. - 28.));
}

} // namespace cpppetsc
} // namespace ae108