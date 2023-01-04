// © 2020 ETH Zurich, Mechanics and Materials Lab
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