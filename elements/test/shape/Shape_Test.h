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

#pragma once

#include "ae108/elements/shape/automatic_gradients.h"
#include "ae108/elements/shape/compute_gradients.h"
#include "ae108/elements/shape/compute_values.h"
#include "ae108/elements/shape/get_points.h"
#include "ae108/elements/tensor/Tensor.h"
#include "ae108/elements/tensor/as_matrix_of_rows.h"
#include <gmock/gmock.h>

namespace ae108 {
namespace elements {
namespace shape {
namespace {

template <class Shape_> struct Shape_Test : ::testing::Test {
  using Shape = Shape_;
};

TYPED_TEST_CASE_P(Shape_Test);

TYPED_TEST_P(Shape_Test, shape_functions_are_correct_at_points) {
  using Shape = typename TestFixture::Shape;

  tensor::Tensor<typename Shape::value_type, Shape::size(), Shape::size()>
      results;
  const auto &points = get_points<Shape>();

  auto point_iterator = points.begin();
  auto result_iterator = results.begin();
  while (point_iterator != points.end() && result_iterator != results.end()) {
    *(result_iterator++) = compute_values<Shape>(*(point_iterator++));
  }

  const auto result_matrix = tensor::as_matrix_of_rows(&results);
  EXPECT_TRUE(result_matrix.Identity().isApprox(result_matrix));
}

TYPED_TEST_P(Shape_Test, gradients_are_correct_at_points) {
  using Shape = typename TestFixture::Shape;

  const auto &points = get_points<Shape>();

  for (const auto &point : points) {
    const auto gradients = compute_gradients<Shape>(point);
    const auto approximation = automatic_gradients<Shape>(point);

    EXPECT_TRUE(tensor::as_matrix_of_rows(&gradients)
                    .isApprox(tensor::as_matrix_of_rows(&approximation)));
  }
}

REGISTER_TYPED_TEST_CASE_P(Shape_Test, shape_functions_are_correct_at_points,
                           gradients_are_correct_at_points);

} // namespace
} // namespace shape
} // namespace elements
} // namespace ae108
