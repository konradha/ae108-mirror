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
#include "ae108/cpppetsc/createTransformInput.h"
#include <gmock/gmock.h>

using testing::SizeIs;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace cpppetsc {

template <class Policy> struct createTransformInput_Test : Test {
  using matrix_type = Matrix<Policy>;
  using size_type = typename Vector<Policy>::size_type;
};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;

TYPED_TEST_CASE(createTransformInput_Test, Policies);

TYPED_TEST(createTransformInput_Test, vector_has_correct_size) {
  using size_type = typename TestFixture::size_type;
  using matrix_type = typename TestFixture::matrix_type;

  const auto rows = size_type{3};
  const auto cols = size_type{7};
  const auto matrix = matrix_type(rows, cols);

  const auto result = createTransformInput(matrix);
  EXPECT_THAT(result.unwrap(), SizeIs(cols));
}

} // namespace cpppetsc
} // namespace ae108