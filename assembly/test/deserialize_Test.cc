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

#include "ae108/assembly/utilities/deserialize.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cppptest/Matchers.h"
#include <Eigen/Core>
#include <array>
#include <gmock/gmock.h>
#include <vector>

using ae108::cppptest::ValueAlmostEq;
using testing::DoubleEq;
using testing::ElementsAre;
using testing::Eq;
using testing::Test;

namespace ae108 {
namespace assembly {
namespace utilities {
namespace {

struct deserialize_Test : Test {};

TEST_F(deserialize_Test, deserializing_vector_works) {
  std::vector<double> values(3);
  const std::vector<double> buffer = {1, 2, 3, 4};

  const auto result = deserialize(buffer.begin(), &values);

  EXPECT_THAT(result, Eq(buffer.begin() + 3));

  EXPECT_THAT(values, ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(3.)));
}

TEST_F(deserialize_Test, deserializing_array_works) {
  auto values = std::array<double, 3>();
  const std::vector<double> buffer = {1, 2, 3, 4};

  const auto result = deserialize(buffer.begin(), &values);

  EXPECT_THAT(result, Eq(buffer.begin() + 3));

  EXPECT_THAT(values, ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(3.)));
}

TEST_F(deserialize_Test, deserializing_eigen_matrix_works) {
  Eigen::Matrix<double, 3, 2> values;
  const std::vector<double> buffer = {1, 2, 3, 4, 5, 6, 7, 8};

  const auto result = deserialize(buffer.begin(), &values);

  EXPECT_THAT(result, Eq(buffer.begin() + 6));

  EXPECT_THAT(values(0, 0), DoubleEq(1.));
  EXPECT_THAT(values(0, 1), DoubleEq(2.));
  EXPECT_THAT(values(1, 0), DoubleEq(3.));
  EXPECT_THAT(values(1, 1), DoubleEq(4.));
  EXPECT_THAT(values(2, 0), DoubleEq(5.));
  EXPECT_THAT(values(2, 1), DoubleEq(6.));
}

TEST_F(deserialize_Test, deserializing_eigen_vector_works) {
  Eigen::Vector3d values;
  const std::vector<double> buffer = {1, 2, 3, 4, 5, 6, 7, 8};

  const auto result = deserialize(buffer.begin(), &values);

  EXPECT_THAT(result, Eq(buffer.begin() + 3));

  EXPECT_THAT(values(0), DoubleEq(1.));
  EXPECT_THAT(values(1), DoubleEq(2.));
  EXPECT_THAT(values(2), DoubleEq(3.));
}

TEST_F(deserialize_Test, deserializing_range_works) {
  std::vector<std::vector<double>> values(2, std::vector<double>(2));
  const std::vector<double> buffer = {1, 2, 3, 4};

  const auto result =
      deserializeRange(buffer.begin(), buffer.end(), values.begin());

  EXPECT_THAT(result, Eq(buffer.begin() + 4));

  EXPECT_THAT(values, ElementsAre(ElementsAre(DoubleEq(1.), DoubleEq(2.)),
                                  ElementsAre(DoubleEq(3.), DoubleEq(4.))));
}

TEST_F(deserialize_Test,
       deserializing_cpppetsc_vector_into_eigen_vectors_works) {
  using vector_type = cpppetsc::Vector<cpppetsc::SequentialComputePolicy>;
  std::vector<Eigen::Matrix<vector_type::value_type, 3, 1>> values(2);
  auto buffer = vector_type::fromList({1., 2., 3., 4., 5., 6.});

  const auto range = buffer.localValues();
  const auto result =
      deserializeRange(range.begin(), range.end(), values.begin());

  EXPECT_THAT(result, Eq(range.end()));

  EXPECT_THAT(values.at(0)(0), ValueAlmostEq(1.));
  EXPECT_THAT(values.at(0)(1), ValueAlmostEq(2.));
  EXPECT_THAT(values.at(0)(2), ValueAlmostEq(3.));
  EXPECT_THAT(values.at(1)(0), ValueAlmostEq(4.));
  EXPECT_THAT(values.at(1)(1), ValueAlmostEq(5.));
  EXPECT_THAT(values.at(1)(2), ValueAlmostEq(6.));
}
} // namespace
} // namespace utilities
} // namespace assembly
} // namespace ae108
