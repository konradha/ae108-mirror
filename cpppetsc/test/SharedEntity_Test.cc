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
#include "ae108/cpppetsc/SharedEntity.h"
#include <gmock/gmock.h>
#include <petscvec.h>

using testing::Eq;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace cpppetsc {
namespace {

template <class Policy> struct SharedEntity_Test : Test {
  /**
   * @brief Create an example vector.
   */
  SharedEntity<Vec, Policy> create() {
    auto vec = Vec{};
    const auto globalSize = PetscInt{10};
    Policy::handleError(
        VecCreateMPI(Policy::communicator(), PETSC_DECIDE, globalSize, &vec));
    return makeSharedEntity<Policy>(vec);
  }
};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;

TYPED_TEST_CASE(SharedEntity_Test, Policies);

TYPED_TEST(SharedEntity_Test, initial_count_is_1) {
  EXPECT_THAT(this->create().count(), Eq(1));
}

TYPED_TEST(SharedEntity_Test, copying_increases_count) {
  const auto value = this->create();
  const SharedEntity<Vec, TypeParam> copy(value);

  EXPECT_THAT(value.count(), Eq(2));
  EXPECT_THAT(copy.count(), Eq(2));
}

TYPED_TEST(SharedEntity_Test, copy_points_to_same_entity) {
  const auto value = this->create();
  const auto copy = value;

  EXPECT_THAT(copy.get(), Eq(value.get()));
}

TYPED_TEST(SharedEntity_Test, assignment_increases_count) {
  const auto value = this->create();
  auto copy = this->create();
  copy = value;

  EXPECT_THAT(value.count(), Eq(2));
  EXPECT_THAT(copy.count(), Eq(2));
}

TYPED_TEST(SharedEntity_Test, assigned_points_to_same_entity) {
  const auto value = this->create();
  auto copy = this->create();
  copy = value;

  EXPECT_THAT(copy.get(), Eq(value.get()));
}

TYPED_TEST(SharedEntity_Test, moving_does_not_increase_count) {
  auto value = this->create();
  const SharedEntity<Vec, TypeParam> copy(std::move(value));

  EXPECT_THAT(copy.count(), Eq(1));
}

TYPED_TEST(SharedEntity_Test, move_assignment_does_not_increase_count) {
  auto value = this->create();
  auto copy = this->create();
  copy = std::move(value);

  EXPECT_THAT(copy.count(), Eq(1));
}

TYPED_TEST(SharedEntity_Test, valid_contents_yield_true) {
  const auto value = this->create();

  EXPECT_THAT(static_cast<bool>(value), Eq(true));
}

TYPED_TEST(SharedEntity_Test, null_contents_yield_false) {
  const auto value = SharedEntity<Vec, TypeParam>(nullptr);

  EXPECT_THAT(static_cast<bool>(value), Eq(false));
}

TYPED_TEST(SharedEntity_Test, destruction_reduces_count) {
  const auto value = this->create();

  {
    const auto copy = value;
    ASSERT_THAT(copy.count(), Eq(2));
  }

  EXPECT_THAT(value.count(), Eq(1));
}

} // namespace
} // namespace cpppetsc
} // namespace ae108