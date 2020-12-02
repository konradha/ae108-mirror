// © 2020 ETH Zurich, Mechanics and Materials Lab
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
#include "ae108/cpppetsc/getName.h"
#include "ae108/cpppetsc/setName.h"
#include <gmock/gmock.h>
#include <string>

using testing::StrEq;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace cpppetsc {
namespace {

template <class Policy> struct setName_Test : Test {};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;

TYPED_TEST_CASE(setName_Test, Policies);

TYPED_TEST(setName_Test, name_of_distributed_vector_can_be_set_and_queried) {
  auto vector = distributed<Vector<TypeParam>>(1);

  const auto name = std::string("test_name");
  setName(name.c_str(), &vector);

  EXPECT_THAT(getName(vector), StrEq(name));
}

} // namespace
} // namespace cpppetsc
} // namespace ae108