// © 2020 ETH Zurich, Mechanics and Materials Lab
// © 2020 California Institute of Technology
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
#include <gmock/gmock.h>

using ae108::cppptest::AlmostEqIfLocal;
using testing::DoubleEq;
using testing::Eq;
using testing::Ge;
using testing::Le;
using testing::Not;
using testing::Pair;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace cpppetsc {
namespace {

template <class Policy> struct Vector_Test : Test {
  using vector_type = Vector<Policy>;
};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;
TYPED_TEST_CASE(Vector_Test, Policies);

TYPED_TEST(Vector_Test, create_vector_of_length_1) {
  typename TestFixture::vector_type vec(1);
  EXPECT_THAT(vec.size(), Eq(1));
}

TYPED_TEST(Vector_Test, local_range_works) {
  typename TestFixture::vector_type vec(3);

  const auto range = vec.localRowRange();
  EXPECT_THAT(range.first, Ge(0));
  EXPECT_THAT(range.second, Ge(range.first));
  EXPECT_THAT(range.second, Le(4));
}

TYPED_TEST(Vector_Test, local_begin_and_end_have_correct_distance) {
  const auto vec = TestFixture::vector_type::fromList({1., 2., 3.});

  const auto range = vec.localRowRange();
  EXPECT_THAT(std::distance(vec.localBegin(), vec.localEnd()),
              Eq(range.second - range.first));
}

TYPED_TEST(Vector_Test, local_begin_and_end_provide_access) {
  const auto vec = TestFixture::vector_type::fromList({1., 2., 3.});

  auto range = vec.localRowRange();
  auto begin = vec.localBegin();

  while (range.first != range.second) {
    ASSERT_THAT(begin, Not(Eq(vec.localEnd())));
    EXPECT_THAT(*begin, Eq(vec(range.first)));
    range.first++;
    begin++;
  }
  EXPECT_THAT(begin, Eq(vec.localEnd()));
}

TYPED_TEST(Vector_Test, range_returns_begin_and_end) {
  const auto vec = TestFixture::vector_type::fromList({1., 2., 3.});

  const auto range = vec.localRange();
  EXPECT_THAT(range.begin(), Eq(vec.localBegin()));
  EXPECT_THAT(range.end(), Eq(vec.localEnd()));
}

TYPED_TEST(Vector_Test, accessing_using_square_brackets_works) {
  const auto vec = TestFixture::vector_type::fromList({1., 2., 3.});

  const auto range = vec.localRowRange();
  for (auto i = range.first; i != range.second; ++i) {
    EXPECT_THAT(vec(i), DoubleEq(vec[i]));
  }
}

TYPED_TEST(Vector_Test, constructing_sequential_copies_works) {
  const auto value = .7;
  distributed<typename TestFixture::vector_type> distributedVector(2, value);
  const auto local =
      TestFixture::vector_type::fromDistributed(distributedVector);

  const auto range = local.unwrap().localRowRange();
  EXPECT_THAT(range, Pair(Eq(0), Eq(2)));
  EXPECT_THAT(local(0), DoubleEq(value));
  EXPECT_THAT(local(1), DoubleEq(value));
}

TYPED_TEST(Vector_Test, move_construction_works) {
  typename TestFixture::vector_type vec(1, 7.);
  auto copy(std::move(vec));

  EXPECT_THAT(copy.size(), Eq(1));
  EXPECT_THAT(copy, AlmostEqIfLocal(0, 7.));
}

TYPED_TEST(Vector_Test, move_assignment_works) {
  typename TestFixture::vector_type vec(1, 7.);
  typename TestFixture::vector_type copy(2, 3.);

  copy = std::move(vec);

  EXPECT_THAT(copy.size(), Eq(1));
  EXPECT_THAT(copy, AlmostEqIfLocal(0, 7.));
}

TYPED_TEST(Vector_Test, swapping_works) {
  typename TestFixture::vector_type vec(1, 7.);
  typename TestFixture::vector_type copy(1, 8.);

  std::swap(vec, copy);

  EXPECT_THAT(copy.size(), Eq(1));
  EXPECT_THAT(copy, AlmostEqIfLocal(0, 7.));
}

TYPED_TEST(Vector_Test, create_vector_of_length_2) {
  typename TestFixture::vector_type vec(2);
  EXPECT_THAT(vec.size(), Eq(2));
}

TYPED_TEST(Vector_Test, from_initializer_list) {
  const auto vec = TestFixture::vector_type::fromList({1, 2, 3});

  EXPECT_THAT(vec.size(), Eq(3));
  EXPECT_THAT(vec, AlmostEqIfLocal(0, 1.));
  EXPECT_THAT(vec, AlmostEqIfLocal(1, 2.));
  EXPECT_THAT(vec, AlmostEqIfLocal(2, 3.));
}

TYPED_TEST(Vector_Test, create_vector_7s_of_length_2) {
  typename TestFixture::vector_type vec(2, 7.);

  EXPECT_THAT(vec, AlmostEqIfLocal(0, 7.));
  EXPECT_THAT(vec, AlmostEqIfLocal(1, 7.));
}

TYPED_TEST(Vector_Test, value_addition_works) {
  typename TestFixture::vector_type vec(2, 3.);

  {
    const auto adder = vec.add();
    if (TypeParam::isPrimaryRank()) {
      adder.element(0, 1.).element(1, 2.);
    }
  }

  EXPECT_THAT(vec, AlmostEqIfLocal(0, 4.));
  EXPECT_THAT(vec, AlmostEqIfLocal(1, 5.));
}

TYPED_TEST(Vector_Test, value_addition_works_in_parallel) {
  typename TestFixture::vector_type vec(2, 3.);

  vec.add().element(0, 1.);

  int size = 0;
  ASSERT_THAT(MPI_Comm_size(TypeParam::communicator(), &size), Eq(0));
  EXPECT_THAT(vec, AlmostEqIfLocal(0, 3. + size));
}

TYPED_TEST(Vector_Test, insertion_works) {
  typename TestFixture::vector_type vec(2);

  vec.replace().element(0, 3.).element(1, 4.);

  EXPECT_THAT(vec, AlmostEqIfLocal(0, 3.));
  EXPECT_THAT(vec, AlmostEqIfLocal(1, 4.));
}

TYPED_TEST(Vector_Test, inserter_offers_size) {
  typename TestFixture::vector_type vec(2);

  EXPECT_THAT(vec.replace().size(), Eq(vec.size()));
}

TYPED_TEST(Vector_Test, insertion_proxy_works) {
  typename TestFixture::vector_type vec(2, 7.);

  vec.replace()(1) = 3.;

  EXPECT_THAT(vec, AlmostEqIfLocal(0, 7.));
  EXPECT_THAT(vec, AlmostEqIfLocal(1, 3.));
}

TYPED_TEST(Vector_Test, addition_proxy_works) {
  typename TestFixture::vector_type vec(2, 7.);

  {
    const auto adder = vec.add();
    if (TypeParam::isPrimaryRank()) {
      adder(1) += 3.;
    }
  }

  EXPECT_THAT(vec, AlmostEqIfLocal(0, 7.));
  EXPECT_THAT(vec, AlmostEqIfLocal(1, 10.));
}

TYPED_TEST(Vector_Test, insertion_overwrites) {
  typename TestFixture::vector_type vec(2, 7.);

  vec.replace().element(0, 5.).element(1, 6.);

  EXPECT_THAT(vec, AlmostEqIfLocal(0, 5.));
  EXPECT_THAT(vec, AlmostEqIfLocal(1, 6.));
}

TYPED_TEST(Vector_Test, wrapping_vectors_works) {
  auto vec = Vec();
  TypeParam::handleError(
      VecCreateMPI(TypeParam::communicator(), PETSC_DECIDE, 1, &vec));

  typename TestFixture::vector_type wrapped_vec(
      UniqueEntity<Vec>(vec, [](Vec) {}));
  wrapped_vec.replace().element(0, .77);

  EXPECT_THAT(wrapped_vec, AlmostEqIfLocal(0, .77));
  VecDestroy(&vec);
}

TYPED_TEST(Vector_Test, filling_vector_works) {
  typename TestFixture::vector_type vec(2, 7.);

  vec.fill(3.);

  EXPECT_THAT(vec, AlmostEqIfLocal(0, 3.));
  EXPECT_THAT(vec, AlmostEqIfLocal(1, 3.));
}

TYPED_TEST(Vector_Test, filling_vector_works_in_inserter) {
  typename TestFixture::vector_type vec(2, 7.);

  vec.replace().fill(3.);

  EXPECT_THAT(vec, AlmostEqIfLocal(0, 3.));
  EXPECT_THAT(vec, AlmostEqIfLocal(1, 3.));
}

TYPED_TEST(Vector_Test, setting_values_to_zero_works) {
  typename TestFixture::vector_type vec(2, 7.);

  vec.setZero();

  EXPECT_THAT(vec, AlmostEqIfLocal(0, 0.));
  EXPECT_THAT(vec, AlmostEqIfLocal(1, 0.));
}

TYPED_TEST(Vector_Test, setting_values_to_zero_works_in_inserter) {
  typename TestFixture::vector_type vec(2, 7.);

  vec.replace().setZero();

  EXPECT_THAT(vec, AlmostEqIfLocal(0, 0.));
  EXPECT_THAT(vec, AlmostEqIfLocal(1, 0.));
}

TYPED_TEST(Vector_Test, duplicating_layout_works) {
  typename TestFixture::vector_type vec(2, 7.);

  const auto newVector = TestFixture::vector_type::fromLayoutOf(vec);

  ASSERT_THAT(newVector.localRowRange(), Eq(vec.localRowRange()));
  EXPECT_THAT(newVector, AlmostEqIfLocal(0, 0.));
  EXPECT_THAT(newVector, AlmostEqIfLocal(1, 0.));
}

TYPED_TEST(Vector_Test, scaling_vector_works) {
  auto vec = TestFixture::vector_type::fromList({7., 8.});

  const auto factor = 2.;
  vec.scale(factor);

  EXPECT_THAT(vec, AlmostEqIfLocal(0, factor * 7.));
  EXPECT_THAT(vec, AlmostEqIfLocal(1, factor * 8.));
}

TYPED_TEST(Vector_Test, adding_alpha_x_works) {
  auto vec_1 = TestFixture::vector_type::fromList({7., 8.});
  const auto vec_2 = TestFixture::vector_type::fromList({2., 3.});

  vec_1.timesAlphaPlusBetaX(2., 3., vec_2);

  ASSERT_THAT(vec_1.size(), Eq(vec_2.size()));
  EXPECT_THAT(vec_1, AlmostEqIfLocal(0, 20.));
  EXPECT_THAT(vec_1, AlmostEqIfLocal(1, 25.));
}

TYPED_TEST(Vector_Test, adding_alpha_x_plus_beta_y_works) {
  auto vec_1 = TestFixture::vector_type::fromList({7., 8.});
  const auto vec_2 = TestFixture::vector_type::fromList({2., 3.});
  const auto vec_3 = TestFixture::vector_type::fromList({-1., -3.});

  vec_1.timesAlphaPlusBetaXPlusGammaY(2., 3., vec_2, 4., vec_3);

  ASSERT_THAT(vec_1.size(), Eq(vec_2.size()));
  EXPECT_THAT(vec_1, AlmostEqIfLocal(0, 16.));
  EXPECT_THAT(vec_1, AlmostEqIfLocal(1, 13.));
}

TYPED_TEST(Vector_Test, adding_Ax_works) {
  auto vec_1 = TestFixture::vector_type::fromList({7., 8.});
  auto vec_2 = TestFixture::vector_type::fromList({2., 3.});
  const auto mat =
      TestFixture::vector_type::matrix_type::fromList({{2., 0.}, {0., 3.}});

  vec_1.addAx(mat, tag<DistributedTag>(std::move(vec_2)));

  ASSERT_THAT(vec_1.size(), Eq(2));
  EXPECT_THAT(vec_1, AlmostEqIfLocal(0, 11.));
  EXPECT_THAT(vec_1, AlmostEqIfLocal(1, 17.));
}

TYPED_TEST(Vector_Test, computing_norm_works) {
  const auto vec = TestFixture::vector_type::fromList({3., 4.});

  EXPECT_THAT(vec.norm(), DoubleEq(5.));
}
} // namespace
} // namespace cpppetsc
} // namespace ae108
