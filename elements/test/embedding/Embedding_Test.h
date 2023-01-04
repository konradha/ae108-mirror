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

#include "ae108/elements/embedding/automatic_jacobian.h"
#include "ae108/elements/embedding/compute_jacobian.h"
#include "ae108/elements/tensor/as_matrix_of_rows.h"
#include <gmock/gmock.h>

namespace ae108 {
namespace elements {
namespace embedding {
namespace {

template <class TestConfiguration> struct Embedding_Test : ::testing::Test {
  using Embedding = typename TestConfiguration::Embedding;
  Embedding embedding = TestConfiguration::create_embedding();
};

TYPED_TEST_CASE_P(Embedding_Test);

TYPED_TEST_P(Embedding_Test, jacobian_is_correct) {
  const auto jacobian = compute_jacobian(
      this->embedding, typename TestFixture::Embedding::ReferencePoint());
  const auto approximation = automatic_jacobian(
      this->embedding, typename TestFixture::Embedding::ReferencePoint());

  EXPECT_TRUE(tensor::as_matrix_of_rows(&jacobian).isApprox(
      tensor::as_matrix_of_rows(&approximation)));
}

REGISTER_TYPED_TEST_CASE_P(Embedding_Test, jacobian_is_correct);

} // namespace
} // namespace embedding
} // namespace elements
} // namespace ae108