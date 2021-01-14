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