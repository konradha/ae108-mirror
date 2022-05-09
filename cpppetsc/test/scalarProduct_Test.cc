// © 2022 ETH Zurich, Mechanics and Materials Lab
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