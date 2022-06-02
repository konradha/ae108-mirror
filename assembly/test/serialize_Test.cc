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

#include "ae108/assembly/utilities/serialize.h"
#include <Eigen/Core>
#include <array>
#include <gmock/gmock.h>
#include <vector>

using testing::DoubleEq;
using testing::ElementsAre;
using testing::Eq;
using testing::Test;

namespace ae108 {
namespace assembly {
namespace utilities {
namespace {

struct serialize_Test : Test {};

TEST_F(serialize_Test, serializing_vector_works) {
  const std::vector<double> values = {1, 2, 3};
  std::vector<double> buffer(4);

  const auto result = serialize(values, buffer.begin());

  EXPECT_THAT(result, Eq(buffer.begin() + 3));

  EXPECT_THAT(buffer, ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(3.),
                                  DoubleEq(0.)));
}

TEST_F(serialize_Test, serializing_array_works) {
  const std::array<double, 3> values = {{1., 2., 3.}};
  std::vector<double> buffer(4);

  const auto result = serialize(values, buffer.begin());

  EXPECT_THAT(result, Eq(buffer.begin() + 3));

  EXPECT_THAT(buffer, ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(3.),
                                  DoubleEq(0.)));
}

TEST_F(serialize_Test, serializing_eigen_matrix_works) {
  Eigen::Matrix<double, 3, 2> values;
  values << 1, 2, 3, 4, 5, 6;
  std::vector<double> buffer(7);

  const auto result = serialize(values, buffer.begin());

  EXPECT_THAT(result, Eq(buffer.begin() + 6));

  EXPECT_THAT(buffer, ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(3.),
                                  DoubleEq(4.), DoubleEq(5.), DoubleEq(6.),
                                  DoubleEq(0.)));
}

TEST_F(serialize_Test, serializing_dynamically_sized_eigen_matrix_works) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> values(3, 2);
  values << 1, 2, 3, 4, 5, 6;
  std::vector<double> buffer(7);

  const auto result = serialize(values, buffer.begin());

  EXPECT_THAT(result, Eq(buffer.begin() + 6));

  EXPECT_THAT(buffer, ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(3.),
                                  DoubleEq(4.), DoubleEq(5.), DoubleEq(6.),
                                  DoubleEq(0.)));
}

TEST_F(serialize_Test, serializing_eigen_vector_works) {
  Eigen::Vector3d values;
  values << 1, 2, 3;
  std::vector<double> buffer(7);

  const auto result = serialize(values, buffer.begin());

  EXPECT_THAT(result, Eq(buffer.begin() + 3));

  EXPECT_THAT(buffer, ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(3.),
                                  DoubleEq(0.), DoubleEq(0.), DoubleEq(0.),
                                  DoubleEq(0.)));
}

TEST_F(serialize_Test, serializing_range_works) {
  const std::vector<std::vector<double>> values({{1, 2}, {3, 4}});
  std::vector<double> buffer(5);

  const auto result =
      serializeRange(values.begin(), values.end(), buffer.begin());

  EXPECT_THAT(result, Eq(buffer.begin() + 4));

  EXPECT_THAT(buffer, ElementsAre(DoubleEq(1.), DoubleEq(2.), DoubleEq(3.),
                                  DoubleEq(4.), DoubleEq(0.)));
}
} // namespace
} // namespace utilities
} // namespace assembly
} // namespace ae108
