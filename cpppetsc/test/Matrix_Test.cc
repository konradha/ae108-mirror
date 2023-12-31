// © 2020 ETH Zurich, Mechanics and Materials Lab
// © 2020 California Institute of Technology
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
#include "ae108/cppptest/Matchers.h"
#include <gmock/gmock.h>

using ae108::cppptest::ScalarEqIfLocal;
using testing::DoubleEq;
using testing::EndsWith;
using testing::Eq;
using testing::Ge;
using testing::HasSubstr;
using testing::Le;
using testing::Pair;
using testing::StartsWith;
using testing::StrEq;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace cpppetsc {
namespace {

template <class Policy> struct Matrix_Test : Test {
  using matrix_type = Matrix<Policy>;
};

using Policies = Types<SequentialComputePolicy, ParallelComputePolicy>;
TYPED_TEST_CASE(Matrix_Test, Policies);

TYPED_TEST(Matrix_Test, size_is_correct) {
  typename TestFixture::matrix_type mat(2, 3);
  EXPECT_THAT(mat.size().first, Eq(2));
  EXPECT_THAT(mat.size().second, Eq(3));
}

TYPED_TEST(Matrix_Test, local_size_is_correct) {
  using matrix_type = typename TestFixture::matrix_type;

  const matrix_type mat(typename matrix_type::LocalRows{2},
                        typename matrix_type::LocalCols{3},
                        typename matrix_type::GlobalRows{PETSC_DETERMINE},
                        typename matrix_type::GlobalCols{PETSC_DETERMINE});

  EXPECT_THAT(mat.localSize(), Pair(2, 3));
}

TYPED_TEST(Matrix_Test, global_size_is_greater_eq_than_local_size) {
  using matrix_type = typename TestFixture::matrix_type;

  const matrix_type mat(typename matrix_type::LocalRows{2},
                        typename matrix_type::LocalCols{3},
                        typename matrix_type::GlobalRows{PETSC_DETERMINE},
                        typename matrix_type::GlobalCols{PETSC_DETERMINE});

  EXPECT_THAT(mat.size(), Pair(Ge(2), Ge(3)));
}

TYPED_TEST(Matrix_Test, default_block_size_is_1) {
  typename TestFixture::matrix_type mat(2, 3);
  EXPECT_THAT(mat.blockSize(), Eq(1));
}

TYPED_TEST(Matrix_Test, move_construction_works) {
  typename TestFixture::matrix_type mat(2, 3);
  auto copy(std::move(mat));
  EXPECT_THAT(copy.size().first, Eq(2));
  EXPECT_THAT(copy.size().second, Eq(3));
}

TYPED_TEST(Matrix_Test, move_assignment_works) {
  typename TestFixture::matrix_type mat(2, 3);
  typename TestFixture::matrix_type copy(1, 2);

  copy = std::move(mat);
  EXPECT_THAT(copy.size().first, Eq(2));
  EXPECT_THAT(copy.size().second, Eq(3));
}

TYPED_TEST(Matrix_Test, local_row_range_works) {
  typename TestFixture::matrix_type mat(2, 3);

  const auto result = mat.localRowRange();
  EXPECT_THAT(result.first, Ge(0));
  EXPECT_THAT(result.second, Ge(result.first));
  EXPECT_THAT(result.second, Le(2));
}

TYPED_TEST(Matrix_Test, local_column_range_works) {
  typename TestFixture::matrix_type mat(2, 3);

  const auto result = mat.localColumnRange();
  EXPECT_THAT(result.first, Ge(0));
  EXPECT_THAT(result.second, Ge(result.first));
  EXPECT_THAT(result.second, Le(3));
}

TYPED_TEST(Matrix_Test, list_initialization_works) {
  const auto mat =
      TestFixture::matrix_type::fromList({{7., 8., 9.}, {10., 11., 12.}});

  ASSERT_THAT(mat.size(), Pair(2, 3));

  EXPECT_THAT(mat, ScalarEqIfLocal(0, 0, 7.));
  EXPECT_THAT(mat, ScalarEqIfLocal(0, 1, 8.));
  EXPECT_THAT(mat, ScalarEqIfLocal(0, 2, 9.));
  EXPECT_THAT(mat, ScalarEqIfLocal(1, 0, 10.));
  EXPECT_THAT(mat, ScalarEqIfLocal(1, 1, 11.));
  EXPECT_THAT(mat, ScalarEqIfLocal(1, 2, 12.));
}

TYPED_TEST(Matrix_Test, list_initialization_works_for_irregular_input) {
  EXPECT_THROW(TestFixture::matrix_type::fromList({{7., 8., 9.}, {10., 11.}}),
               InvalidParametersException);
}

TYPED_TEST(Matrix_Test, list_initialization_throws_for_empty_input) {
  EXPECT_THROW(TestFixture::matrix_type::fromList({{7., 8., 9.}, {}}),
               InvalidParametersException);
}

TYPED_TEST(Matrix_Test, adding_single_element_works) {
  auto mat =
      TestFixture::matrix_type::fromList({{7., 8., 9.}, {10., 11., 12.}});

  {
    const auto adder = mat.assemblyView().add();
    if (TypeParam::isPrimaryRank()) {
      adder.element(1, 2, 7.);
    }
  }

  EXPECT_THAT(mat, ScalarEqIfLocal(1, 2, 12. + 7.));
}

TYPED_TEST(Matrix_Test, setting_values_works) {
  typename TestFixture::matrix_type mat(2, 3);

  {
    const auto replacer = mat.assemblyView().replace();
    if (TypeParam::isPrimaryRank()) {
      replacer.elements({1}, {0, 2}, {11., 12.});
    }
  }

  EXPECT_THAT(mat, ScalarEqIfLocal(0, 0, 0.));
  EXPECT_THAT(mat, ScalarEqIfLocal(0, 1, 0.));
  EXPECT_THAT(mat, ScalarEqIfLocal(0, 2, 0.));
  EXPECT_THAT(mat, ScalarEqIfLocal(1, 0, 11.));
  EXPECT_THAT(mat, ScalarEqIfLocal(1, 1, 0.));
  EXPECT_THAT(mat, ScalarEqIfLocal(1, 2, 12.));
}

TYPED_TEST(Matrix_Test, inserter_offers_rows) {
  typename TestFixture::matrix_type mat(2, 3);
  EXPECT_THAT(mat.assemblyView().add().rows(), 2);
}

TYPED_TEST(Matrix_Test, inserter_offers_cols) {
  typename TestFixture::matrix_type mat(2, 3);
  EXPECT_THAT(mat.assemblyView().add().cols(), 3);
}

TYPED_TEST(Matrix_Test, setting_values_via_proxy_works) {
  typename TestFixture::matrix_type mat(1, 2);

  {
    const auto replacer = mat.assemblyView().replace();
    if (TypeParam::isPrimaryRank()) {
      replacer(0, 1) = 7.;
    }
  }

  EXPECT_THAT(mat, ScalarEqIfLocal(0, 0, 0.));
  EXPECT_THAT(mat, ScalarEqIfLocal(0, 1, 7.));
}

TYPED_TEST(Matrix_Test, adding_values_via_proxy_works) {
  auto mat = TestFixture::matrix_type::fromList({{7., 7.}});

  {
    const auto adder = mat.assemblyView().add();
    if (TypeParam::isPrimaryRank()) {
      adder(0, 1) += 3.;
    }
  }

  EXPECT_THAT(mat, ScalarEqIfLocal(0, 0, 7.));
  EXPECT_THAT(mat, ScalarEqIfLocal(0, 1, 10.));
}

TYPED_TEST(Matrix_Test, setting_values_via_allocating_view_works) {
  typename TestFixture::matrix_type mat(2, 2);

  {
    const auto replacer = mat.preallocatedAssemblyView(1).replace();
    if (TypeParam::isPrimaryRank()) {
      replacer(0, 0) = 7.;
      replacer(1, 1) = 8.;
    }
  }

  EXPECT_THAT(mat, ScalarEqIfLocal(0, 0, 7.));
  EXPECT_THAT(mat, ScalarEqIfLocal(0, 1, 0.));
  EXPECT_THAT(mat, ScalarEqIfLocal(1, 0, 0.));
  EXPECT_THAT(mat, ScalarEqIfLocal(1, 1, 8.));
}

TYPED_TEST(Matrix_Test, wrapping_matrices_works) {
  Mat mat;

  TypeParam::handleError(MatCreate(TypeParam::communicator(), &mat));
  TypeParam::handleError(
      MatSetSizes(mat, PETSC_DETERMINE, PETSC_DETERMINE, 1, 1));
  TypeParam::handleError(MatSetUp(mat));

  typename TestFixture::matrix_type wrapped_mat(
      UniqueEntity<Mat>(mat, [](Mat) {}));
  wrapped_mat.assemblyView().replace().element(0, 0, .77);

  EXPECT_THAT(wrapped_mat, ScalarEqIfLocal(0, 0, .77));
  MatDestroy(&mat);
}

TYPED_TEST(Matrix_Test, zeroing_lines_works) {
  auto mat = TestFixture::matrix_type::fromList(
      {{1., 2., 3.}, {4., 5., 6.}, {7., 8., 9.}});

  mat.replaceRowsByEye({0, 2});

  EXPECT_THAT(mat, ScalarEqIfLocal(0, 0, 1.));
  EXPECT_THAT(mat, ScalarEqIfLocal(0, 1, 0.));
  EXPECT_THAT(mat, ScalarEqIfLocal(0, 2, 0.));
  EXPECT_THAT(mat, ScalarEqIfLocal(1, 0, 4.));
  EXPECT_THAT(mat, ScalarEqIfLocal(1, 1, 5.));
  EXPECT_THAT(mat, ScalarEqIfLocal(1, 2, 6.));
  EXPECT_THAT(mat, ScalarEqIfLocal(2, 0, 0.));
  EXPECT_THAT(mat, ScalarEqIfLocal(2, 1, 0.));
  EXPECT_THAT(mat, ScalarEqIfLocal(2, 2, 1.));
}

TYPED_TEST(Matrix_Test, zeroing_all_values_works) {
  auto mat = TestFixture::matrix_type::fromList({
      {1., 2.},
      {4., 5.},
  });

  mat.setZero();

  EXPECT_THAT(mat, ScalarEqIfLocal(0, 0, 0.));
  EXPECT_THAT(mat, ScalarEqIfLocal(0, 1, 0.));
  EXPECT_THAT(mat, ScalarEqIfLocal(1, 0, 0.));
  EXPECT_THAT(mat, ScalarEqIfLocal(1, 1, 0.));
}

TYPED_TEST(Matrix_Test, zeroing_all_values_works_for_inserter) {
  auto mat = TestFixture::matrix_type::fromList({
      {1., 2.},
      {4., 5.},
  });

  mat.assemblyView().replace().setZero();

  EXPECT_THAT(mat, ScalarEqIfLocal(0, 0, 0.));
  EXPECT_THAT(mat, ScalarEqIfLocal(0, 1, 0.));
  EXPECT_THAT(mat, ScalarEqIfLocal(1, 0, 0.));
  EXPECT_THAT(mat, ScalarEqIfLocal(1, 1, 0.));
}

TYPED_TEST(Matrix_Test, multadding_matrix_works) {
  auto mat_1 = TestFixture::matrix_type::fromList({
      {1., 2.},
      {4., 5.},
  });

  const auto mat_2 = TestFixture::matrix_type::fromList({
      {1., 3.},
      {7., 9.},
  });

  mat_1.addAlphaX(-2., mat_2);

  EXPECT_THAT(mat_1, ScalarEqIfLocal(0, 0, -1.));
  EXPECT_THAT(mat_1, ScalarEqIfLocal(0, 1, -4.));
  EXPECT_THAT(mat_1, ScalarEqIfLocal(1, 0, -10.));
  EXPECT_THAT(mat_1, ScalarEqIfLocal(1, 1, -13.));
}

TYPED_TEST(Matrix_Test, scaling_matrix_works) {
  auto mat = TestFixture::matrix_type::fromList({
      {1., 2.},
      {4., 5.},
  });

  const auto factor = 2.;
  mat.scale(factor);

  EXPECT_THAT(mat, ScalarEqIfLocal(0, 0, 2.));
  EXPECT_THAT(mat, ScalarEqIfLocal(0, 1, 4.));
  EXPECT_THAT(mat, ScalarEqIfLocal(1, 0, 8.));
  EXPECT_THAT(mat, ScalarEqIfLocal(1, 1, 10.));
}

TYPED_TEST(Matrix_Test, duplicating_layout_works) {
  const auto mat = TestFixture::matrix_type::fromList({
      {1., 2., 3.},
      {4., 5., 6.},
  });

  auto duplicate = TestFixture::matrix_type::fromLayoutOf(mat);

  ASSERT_THAT(duplicate.size().first, Eq(mat.size().first));
  ASSERT_THAT(duplicate.size().second, Eq(mat.size().second));

  EXPECT_THAT(duplicate, ScalarEqIfLocal(0, 0, 0.));
  EXPECT_THAT(duplicate, ScalarEqIfLocal(0, 1, 0.));
  EXPECT_THAT(duplicate, ScalarEqIfLocal(0, 2, 0.));
  EXPECT_THAT(duplicate, ScalarEqIfLocal(1, 0, 0.));
  EXPECT_THAT(duplicate, ScalarEqIfLocal(1, 1, 0.));
  EXPECT_THAT(duplicate, ScalarEqIfLocal(1, 2, 0.));
}

TYPED_TEST(Matrix_Test, cloning_works) {
  const auto mat = TestFixture::matrix_type::fromList({
      {1., 2., 3.},
      {4., 5., 6.},
  });

  auto duplicate = TestFixture::matrix_type::clone(mat);

  ASSERT_THAT(duplicate.size().first, Eq(mat.size().first));
  ASSERT_THAT(duplicate.size().second, Eq(mat.size().second));

  EXPECT_THAT(duplicate, ScalarEqIfLocal(0, 0, 1.));
  EXPECT_THAT(duplicate, ScalarEqIfLocal(0, 1, 2.));
  EXPECT_THAT(duplicate, ScalarEqIfLocal(0, 2, 3.));
  EXPECT_THAT(duplicate, ScalarEqIfLocal(1, 0, 4.));
  EXPECT_THAT(duplicate, ScalarEqIfLocal(1, 1, 5.));
  EXPECT_THAT(duplicate, ScalarEqIfLocal(1, 2, 6.));
}

TYPED_TEST(Matrix_Test, inverting_block_diagonal_works) {
  const auto A = TestFixture::matrix_type::fromList({
      {3., 2.},
      {4., 5.},
  });

  ASSERT_THAT(A.blockSize(), Eq(1));

  const auto result = TestFixture::matrix_type::fromInvertedBlockDiagonalOf(A);

  ASSERT_THAT(result.size().first, Eq(2));
  ASSERT_THAT(result.size().second, Eq(2));

  EXPECT_THAT(result, ScalarEqIfLocal(0, 0, 1. / 3.));
  EXPECT_THAT(result, ScalarEqIfLocal(0, 1, 0.));
  EXPECT_THAT(result, ScalarEqIfLocal(1, 0, 0.));
  EXPECT_THAT(result, ScalarEqIfLocal(1, 1, 1. / 5.));
}

TYPED_TEST(Matrix_Test, double_product_works) {
  const auto A = TestFixture::matrix_type::fromList({{0., 1.}});
  const auto B = TestFixture::matrix_type::fromList({
      {1., 2., 3.},
      {4., 5., 6.},
  });

  const auto result = TestFixture::matrix_type::fromProduct(A, B);

  ASSERT_THAT(result.size().first, Eq(1));
  ASSERT_THAT(result.size().second, Eq(3));

  EXPECT_THAT(result, ScalarEqIfLocal(0, 0, 4.));
  EXPECT_THAT(result, ScalarEqIfLocal(0, 1, 5.));
  EXPECT_THAT(result, ScalarEqIfLocal(0, 2, 6.));
}

TYPED_TEST(Matrix_Test, transpose_matrix_product_works) {
  const auto A = TestFixture::matrix_type::fromList({{0.}, {1.}});
  const auto B = TestFixture::matrix_type::fromList({
      {1., 2., 3.},
      {4., 5., 6.},
  });

  const auto result = TestFixture::matrix_type::fromAtB(A, B);

  ASSERT_THAT(result.size().first, Eq(1));
  ASSERT_THAT(result.size().second, Eq(3));

  EXPECT_THAT(result, ScalarEqIfLocal(0, 0, 4.));
  EXPECT_THAT(result, ScalarEqIfLocal(0, 1, 5.));
  EXPECT_THAT(result, ScalarEqIfLocal(0, 2, 6.));
}

TYPED_TEST(Matrix_Test, triple_product_works) {
  const auto A = TestFixture::matrix_type::fromList({{0., 1.}});
  const auto B = TestFixture::matrix_type::fromList({
      {1., 2., 3.},
      {4., 5., 6.},
  });
  const auto C = TestFixture::matrix_type::fromList({{1.}, {0.}, {0.}});

  const auto result = TestFixture::matrix_type::fromProduct(A, B, C);

  ASSERT_THAT(result.size().first, Eq(1));
  ASSERT_THAT(result.size().second, Eq(1));

  EXPECT_THAT(result, ScalarEqIfLocal(0, 0, 4.));
}

TYPED_TEST(Matrix_Test, PtAP_works) {
  const auto P = TestFixture::matrix_type::fromList({{0.}, {1.}});
  const auto A = TestFixture::matrix_type::fromList({
      {1., 2.},
      {4., 5.},
  });

  const auto result = TestFixture::matrix_type::fromPtAP(P, A);

  ASSERT_THAT(result.size().first, Eq(1));
  ASSERT_THAT(result.size().second, Eq(1));

  EXPECT_THAT(result, ScalarEqIfLocal(0, 0, 5.));
}

TYPED_TEST(Matrix_Test, PAPt_works) {
  const auto P = TestFixture::matrix_type::fromList({{0., 1.}});
  const auto A = TestFixture::matrix_type::fromList({
      {1., 2.},
      {4., 5.},
  });

  const auto result = TestFixture::matrix_type::fromPAPt(P, A);

  ASSERT_THAT(result.size().first, Eq(1));
  ASSERT_THAT(result.size().second, Eq(1));

  EXPECT_THAT(result, ScalarEqIfLocal(0, 0, 5.));
}

TYPED_TEST(Matrix_Test, computing_the_norm_works) {
  const auto mat = TestFixture::matrix_type::fromList({
      {3., 0.},
      {0., 4.},
  });

  EXPECT_THAT(mat.norm(), DoubleEq(5.));
}

/**
 * @brief Converts the parameter to a string using a stringstream and
 * operator<<.
 */
template <class T> std::string toString(const T &t) {
  std::stringstream stream;
  stream << t;
  return stream.str();
}

TYPED_TEST(Matrix_Test, writing_to_stream_uses_square_brackets) {
  const auto mat = TestFixture::matrix_type::fromList({
      {3., 0.},
      {0., 4.},
  });

  const auto result = toString(mat);
  const auto hasRows =
      bool{mat.localRowRange().second - mat.localRowRange().first > 0};

  if (hasRows) {
    EXPECT_THAT(result, StartsWith("[["));
    EXPECT_THAT(result, EndsWith("]]"));
  } else {
    EXPECT_THAT(result, StrEq("[]"));
  }
}

TYPED_TEST(Matrix_Test, local_values_are_written_to_stream) {
  const auto mat = TestFixture::matrix_type::fromList({
      {3., 0.},
      {0., 4.},
  });

  const auto result = toString(mat);

  const auto range = mat.localRowRange();
  const auto size = mat.size();
  for (auto row = range.first; row < range.second; ++row) {
    for (auto col = decltype(size.second){0}; col < size.second; ++col) {
      EXPECT_THAT(result, HasSubstr(toString(mat(row, col))));
    }
  }
}

} // namespace
} // namespace cpppetsc
} // namespace ae108
