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

#include "ae108/elements/mesh/refine_segment_mesh.h"
#include <gmock/gmock.h>

using testing::DoubleEq;
using testing::Eq;
using testing::IsEmpty;
using testing::SizeIs;
using testing::Test;

namespace ae108 {
namespace elements {
namespace mesh {
namespace {

struct refine_segment_mesh_Test : Test {};

TEST_F(refine_segment_mesh_Test,
       returns_correct_connectivity_for_no_refinement_2D) {

  using Point = tensor::Tensor<double, 2>;
  const elements::mesh::Mesh<Point> unrefined_mesh{{{0, 1}},
                                                   {{0., 0.}, {1., 0.}}};

  const auto refined_mesh = refine_segment_mesh(unrefined_mesh, 2);

  ASSERT_THAT(refined_mesh.connectivity(), SizeIs(1));
  EXPECT_THAT(refined_mesh.connectivity().at(0).at(0), Eq(0));
  EXPECT_THAT(refined_mesh.connectivity().at(0).at(1), Eq(1));
}

TEST_F(refine_segment_mesh_Test,
       returns_correct_positions_for_no_refinement_2D) {

  using Point = tensor::Tensor<double, 2>;
  const elements::mesh::Mesh<Point> unrefined_mesh{{{0, 1}},
                                                   {{0., 0.}, {1., 0.}}};

  const auto refined_mesh = refine_segment_mesh(unrefined_mesh, 2);

  ASSERT_THAT(refined_mesh.number_of_positions(), Eq(2));
  EXPECT_THAT(refined_mesh.position_of_vertex(0).at(0), DoubleEq(0.));
  EXPECT_THAT(refined_mesh.position_of_vertex(0).at(1), DoubleEq(0.));
  EXPECT_THAT(refined_mesh.position_of_vertex(1).at(0), DoubleEq(1.));
  EXPECT_THAT(refined_mesh.position_of_vertex(1).at(1), DoubleEq(0.));
}

TEST_F(refine_segment_mesh_Test,
       returns_correct_connectivity_for_no_refinement_3D) {

  using Point = tensor::Tensor<double, 3>;
  const elements::mesh::Mesh<Point> unrefined_mesh{
      {{0, 1}}, {{0., 0., 0.}, {1., 0., 0.}}};

  const auto refined_mesh = refine_segment_mesh(unrefined_mesh, 2);

  ASSERT_THAT(refined_mesh.connectivity(), SizeIs(1));
  EXPECT_THAT(refined_mesh.connectivity().at(0).at(0), Eq(0));
  EXPECT_THAT(refined_mesh.connectivity().at(0).at(1), Eq(1));
}

TEST_F(refine_segment_mesh_Test,
       returns_correct_positions_for_no_refinement_3D) {

  using Point = tensor::Tensor<double, 3>;
  const elements::mesh::Mesh<Point> unrefined_mesh{
      {{0, 1}}, {{0., 0., 0.}, {1., 0., 0.}}};

  const auto refined_mesh = refine_segment_mesh(unrefined_mesh, 2);

  ASSERT_THAT(refined_mesh.number_of_positions(), Eq(2));
  EXPECT_THAT(refined_mesh.position_of_vertex(0).at(0), DoubleEq(0.));
  EXPECT_THAT(refined_mesh.position_of_vertex(0).at(1), DoubleEq(0.));
  EXPECT_THAT(refined_mesh.position_of_vertex(0).at(2), DoubleEq(0.));
  EXPECT_THAT(refined_mesh.position_of_vertex(1).at(0), DoubleEq(1.));
  EXPECT_THAT(refined_mesh.position_of_vertex(1).at(1), DoubleEq(0.));
  EXPECT_THAT(refined_mesh.position_of_vertex(1).at(2), DoubleEq(0.));
}

TEST_F(refine_segment_mesh_Test,
       returns_correct_connectivity_for_single_refinement_2D) {

  using Point = tensor::Tensor<double, 2>;
  const elements::mesh::Mesh<Point> unrefined_mesh{{{0, 1}},
                                                   {{0., 0.}, {1., 0.}}};

  const auto refined_mesh = refine_segment_mesh(unrefined_mesh, 0.6);

  ASSERT_THAT(refined_mesh.connectivity(), SizeIs(2));
  EXPECT_THAT(refined_mesh.connectivity().at(0).at(0), Eq(0));
  EXPECT_THAT(refined_mesh.connectivity().at(0).at(1), Eq(2));
  EXPECT_THAT(refined_mesh.connectivity().at(1).at(0), Eq(2));
  EXPECT_THAT(refined_mesh.connectivity().at(1).at(1), Eq(1));
}

TEST_F(refine_segment_mesh_Test,
       returns_correct_positions_for_single_refinement_2D) {

  using Point = tensor::Tensor<double, 2>;
  const elements::mesh::Mesh<Point> unrefined_mesh{{{0, 1}},
                                                   {{0., 0.}, {1., 0.}}};

  const auto refined_mesh = refine_segment_mesh(unrefined_mesh, 0.6);

  ASSERT_THAT(refined_mesh.number_of_positions(), Eq(3));
  EXPECT_THAT(refined_mesh.position_of_vertex(0).at(0), DoubleEq(0.));
  EXPECT_THAT(refined_mesh.position_of_vertex(0).at(1), DoubleEq(0.));
  EXPECT_THAT(refined_mesh.position_of_vertex(1).at(0), DoubleEq(1.));
  EXPECT_THAT(refined_mesh.position_of_vertex(1).at(1), DoubleEq(0.));
  EXPECT_THAT(refined_mesh.position_of_vertex(2).at(0), DoubleEq(0.5));
  EXPECT_THAT(refined_mesh.position_of_vertex(2).at(1), DoubleEq(0.));
}

TEST_F(refine_segment_mesh_Test,
       returns_correct_connectivity_for_double_refinement_2D) {

  using Point = tensor::Tensor<double, 2>;
  const elements::mesh::Mesh<Point> unrefined_mesh{{{0, 1}},
                                                   {{0., 0.}, {1., 0.}}};

  const auto refined_mesh = refine_segment_mesh(unrefined_mesh, 0.4);

  ASSERT_THAT(refined_mesh.connectivity(), SizeIs(3));
  EXPECT_THAT(refined_mesh.connectivity().at(0).at(0), Eq(0));
  EXPECT_THAT(refined_mesh.connectivity().at(0).at(1), Eq(2));
  EXPECT_THAT(refined_mesh.connectivity().at(1).at(0), Eq(2));
  EXPECT_THAT(refined_mesh.connectivity().at(1).at(1), Eq(3));
  EXPECT_THAT(refined_mesh.connectivity().at(2).at(0), Eq(3));
  EXPECT_THAT(refined_mesh.connectivity().at(2).at(1), Eq(1));
}

TEST_F(refine_segment_mesh_Test,
       returns_correct_positions_for_double_refinement_2D) {

  using Point = tensor::Tensor<double, 2>;
  const elements::mesh::Mesh<Point> unrefined_mesh{{{0, 1}},
                                                   {{0., 0.}, {1., 0.}}};

  const auto refined_mesh = refine_segment_mesh(unrefined_mesh, 0.4);

  ASSERT_THAT(refined_mesh.number_of_positions(), Eq(4));
  EXPECT_THAT(refined_mesh.position_of_vertex(0).at(0), DoubleEq(0.));
  EXPECT_THAT(refined_mesh.position_of_vertex(0).at(1), DoubleEq(0.));
  EXPECT_THAT(refined_mesh.position_of_vertex(1).at(0), DoubleEq(1.));
  EXPECT_THAT(refined_mesh.position_of_vertex(1).at(1), DoubleEq(0.));
  EXPECT_THAT(refined_mesh.position_of_vertex(2).at(0), DoubleEq(1. / 3));
  EXPECT_THAT(refined_mesh.position_of_vertex(2).at(1), DoubleEq(0.));
  EXPECT_THAT(refined_mesh.position_of_vertex(3).at(0), DoubleEq(2. / 3));
  EXPECT_THAT(refined_mesh.position_of_vertex(3).at(1), DoubleEq(0.));
}

TEST_F(refine_segment_mesh_Test,
       returns_correct_connectivity_for_single_refinement_3D) {

  using Point = tensor::Tensor<double, 3>;
  const elements::mesh::Mesh<Point> unrefined_mesh{
      {{0, 1}}, {{0., 0., 0.}, {1., 0., 0.}}};

  const auto refined_mesh = refine_segment_mesh(unrefined_mesh, 0.6);

  ASSERT_THAT(refined_mesh.connectivity(), SizeIs(2));
  EXPECT_THAT(refined_mesh.connectivity().at(0).at(0), Eq(0));
  EXPECT_THAT(refined_mesh.connectivity().at(0).at(1), Eq(2));
  EXPECT_THAT(refined_mesh.connectivity().at(1).at(0), Eq(2));
  EXPECT_THAT(refined_mesh.connectivity().at(1).at(1), Eq(1));
}

TEST_F(refine_segment_mesh_Test,
       returns_correct_positions_for_single_refinement_3D) {

  using Point = tensor::Tensor<double, 3>;
  const elements::mesh::Mesh<Point> unrefined_mesh{
      {{0, 1}}, {{0., 0., 0.}, {1., 0., 0.}}};

  const auto refined_mesh = refine_segment_mesh(unrefined_mesh, 0.6);

  ASSERT_THAT(refined_mesh.number_of_positions(), Eq(3));
  EXPECT_THAT(refined_mesh.position_of_vertex(0).at(0), DoubleEq(0.));
  EXPECT_THAT(refined_mesh.position_of_vertex(0).at(1), DoubleEq(0.));
  EXPECT_THAT(refined_mesh.position_of_vertex(0).at(2), DoubleEq(0.));
  EXPECT_THAT(refined_mesh.position_of_vertex(1).at(0), DoubleEq(1.));
  EXPECT_THAT(refined_mesh.position_of_vertex(1).at(1), DoubleEq(0.));
  EXPECT_THAT(refined_mesh.position_of_vertex(1).at(2), DoubleEq(0.));
  EXPECT_THAT(refined_mesh.position_of_vertex(2).at(0), DoubleEq(0.5));
  EXPECT_THAT(refined_mesh.position_of_vertex(2).at(1), DoubleEq(0.));
  EXPECT_THAT(refined_mesh.position_of_vertex(2).at(2), DoubleEq(0.));
}

TEST_F(refine_segment_mesh_Test,
       returns_correct_number_of_segments_for_complex_refinement_2D) {

  using Point = tensor::Tensor<double, 2>;
  const elements::mesh::Mesh<Point> unrefined_mesh{
      {{0, 1}, {1, 2}, {2, 3}, {3, 0}},
      {{0., 0.}, {1., 0.}, {1., 1.}, {0., 1.}}};

  const auto refined_mesh = refine_segment_mesh(unrefined_mesh, 0.1);

  ASSERT_THAT(refined_mesh.connectivity(), SizeIs(40));
}

TEST_F(refine_segment_mesh_Test,
       returns_correct_number_of_positions_for_complex_refinement_2D) {

  using Point = tensor::Tensor<double, 2>;
  const elements::mesh::Mesh<Point> unrefined_mesh{
      {{0, 1}, {1, 2}, {2, 3}, {3, 0}},
      {{0., 0.}, {1., 0.}, {1., 1.}, {0., 1.}}};

  const auto refined_mesh = refine_segment_mesh(unrefined_mesh, 0.1);

  ASSERT_THAT(refined_mesh.number_of_positions(), Eq(40));
}

} // namespace
} // namespace mesh
} // namespace elements
} // namespace ae108