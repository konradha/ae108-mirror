// © 2021 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/meshing/construct_periodic_point_cloud.h"
#include <gmock/gmock.h>

using testing::DoubleEq;
using testing::ElementsAre;
using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
struct construct_periodic_point_cloud_Test : Test {};

TEST_F(construct_periodic_point_cloud_Test,
       returns_correct_nonsymmetric_3D_point_cloud) {

  const auto point_cloud = construct_periodic_point_cloud<3>(
      {{{1, 0, 0}, {0, 2, 0}, {0, 0, 3}}}, {0, 0, 0}, {1, 0});

  ASSERT_THAT(point_cloud.size(), Eq(8));
  ASSERT_THAT(point_cloud[0],
              ElementsAre(DoubleEq(1), DoubleEq(2), DoubleEq(3)));
  ASSERT_THAT(point_cloud[7],
              ElementsAre(DoubleEq(0), DoubleEq(0), DoubleEq(0)));
};

TEST_F(construct_periodic_point_cloud_Test,
       returns_correct_3D_point_cloud_around_zero_origin) {

  const auto point_cloud = construct_periodic_point_cloud<3>(
      {{{1, 0, 0}, {0, 2, 0}, {0, 0, 3}}}, {0, 0, 0});

  ASSERT_THAT(point_cloud.size(), Eq(27));
  ASSERT_THAT(point_cloud[0],
              ElementsAre(DoubleEq(-1), DoubleEq(-2), DoubleEq(-3)));
  ASSERT_THAT(point_cloud[26],
              ElementsAre(DoubleEq(1), DoubleEq(2), DoubleEq(3)));
};

TEST_F(construct_periodic_point_cloud_Test,
       returns_correct_3D_point_cloud_around_nonzero_origin) {

  const auto point_cloud = construct_periodic_point_cloud<3>(
      {{{1, 0, 0}, {0, 2, 0}, {0, 0, 3}}}, {1, 2, 3});

  ASSERT_THAT(point_cloud.size(), Eq(27));
  ASSERT_THAT(point_cloud[0],
              ElementsAre(DoubleEq(0), DoubleEq(0), DoubleEq(0)));
  ASSERT_THAT(point_cloud[26],
              ElementsAre(DoubleEq(2), DoubleEq(4), DoubleEq(6)));
};

TEST_F(construct_periodic_point_cloud_Test,
       returns_correct_2D_point_cloud_around_zero_origin) {

  const auto point_cloud =
      construct_periodic_point_cloud<2>({{{1, 0}, {0, 2}}}, {0, 0});

  ASSERT_THAT(point_cloud.size(), Eq(9));
  ASSERT_THAT(point_cloud[0], ElementsAre(DoubleEq(-1), DoubleEq(-2)));
  ASSERT_THAT(point_cloud[8], ElementsAre(DoubleEq(1), DoubleEq(2)));
};

TEST_F(construct_periodic_point_cloud_Test,
       returns_correct_1D_point_cloud_around_zero_origin) {

  const auto point_cloud = construct_periodic_point_cloud<1>({{{2}}}, {0});

  ASSERT_THAT(point_cloud.size(), Eq(3));
  ASSERT_THAT(point_cloud[0], ElementsAre(DoubleEq(-2)));
  ASSERT_THAT(point_cloud[2], ElementsAre(DoubleEq(2)));
};

} // namespace meshing
} // namespace ae108