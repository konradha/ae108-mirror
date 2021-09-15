// © 2021 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/cpppetsc/Mesh.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"
#include "ae108/cpppetsc/createTransformInput.h"
#include "ae108/cppptest/Matchers.h"
#include "ae108/solve/InconsistentBoundaryConditionsException.h"
#include "ae108/solve/InvalidVertexException.h"
#include "ae108/solve/boundaryConditionsToTransform.h"
#include <array>
#include <gmock/gmock.h>

using ae108::cppptest::ScalarEqIfLocal;
using testing::DoubleEq;
using testing::Test;

namespace ae108 {
namespace solve {
namespace {

template <class Policy> struct boundaryConditionsToTransform_Test : Test {
  using mesh_type = cpppetsc::Mesh<Policy>;
  using size_type = typename mesh_type::size_type;
  using value_type = typename mesh_type::value_type;

  static constexpr size_type numberOfVertices = size_type{2};
  static constexpr size_type numberOfElements = size_type{2};
  static constexpr size_type dimension = size_type{1};
  static constexpr size_type dofPerVertex = size_type{3};
  static constexpr size_type verticesPerElement = size_type{1};

  using Connectivity =
      std::array<std::array<size_type, verticesPerElement>, numberOfElements>;
  Connectivity connectivity = {{
      {{0}},
      {{1}},
  }}; // two elements with a single vertex each
  mesh_type mesh = mesh_type::fromConnectivity(dimension, connectivity,
                                               numberOfVertices, dofPerVertex);
};

using Policies = testing::Types<cpppetsc::SequentialComputePolicy,
                                cpppetsc::ParallelComputePolicy>;

TYPED_TEST_CASE(boundaryConditionsToTransform_Test, Policies);

TYPED_TEST(boundaryConditionsToTransform_Test,
           returns_identity_if_no_conditions) {
  using size_type = typename TestFixture::size_type;
  const auto result = boundaryConditionsToTransform({}, this->mesh);

  for (size_type i = 0; i < this->numberOfVertices * this->dofPerVertex; ++i)
    for (size_type j = 0; j < this->numberOfVertices * this->dofPerVertex; ++j)
      EXPECT_THAT(result.matrix, ScalarEqIfLocal(i, j, i == j ? 1. : 0.));

  EXPECT_THAT(result.vector.unwrap().norm(), DoubleEq(0.));
}

TYPED_TEST(boundaryConditionsToTransform_Test, dof_can_be_set_to_fixed_value) {
  using size_type = typename TestFixture::size_type;
  using value_type = typename TestFixture::value_type;

  constexpr auto vertex = size_type{1};
  constexpr auto dof = size_type{2};
  constexpr auto value = value_type{.7};

  const auto result =
      boundaryConditionsToTransform({{{vertex, dof}, {}, value}}, this->mesh);
  const auto image = this->mesh.toCanonicalOrder(
      apply(result, createTransformInput(result.matrix)));

  EXPECT_THAT(image.unwrap(),
              ScalarEqIfLocal(vertex * this->dofPerVertex + dof, value));
}

TYPED_TEST(boundaryConditionsToTransform_Test,
           two_dofs_can_be_set_to_fixed_value) {
  using size_type = typename TestFixture::size_type;
  using value_type = typename TestFixture::value_type;

  constexpr auto vertex_1 = size_type{1};
  constexpr auto dof_1 = size_type{2};
  constexpr auto value_1 = value_type{.7};
  constexpr auto vertex_2 = size_type{0};
  constexpr auto dof_2 = size_type{1};
  constexpr auto value_2 = value_type{.3};

  const auto result = boundaryConditionsToTransform(
      {
          {{vertex_1, dof_1}, {}, value_1},
          {{vertex_2, dof_2}, {}, value_2},
      },
      this->mesh);
  const auto image = this->mesh.toCanonicalOrder(
      apply(result, createTransformInput(result.matrix)));

  EXPECT_THAT(image.unwrap(),
              ScalarEqIfLocal(vertex_1 * this->dofPerVertex + dof_1, value_1));
  EXPECT_THAT(image.unwrap(),
              ScalarEqIfLocal(vertex_2 * this->dofPerVertex + dof_2, value_2));
}

TYPED_TEST(boundaryConditionsToTransform_Test,
           dof_can_be_set_to_twice_other_dof) {
  using size_type = typename TestFixture::size_type;
  using value_type = typename TestFixture::value_type;

  constexpr auto targetVertex = size_type{0};
  constexpr auto targetDof = size_type{1};
  constexpr auto sourceVertex = size_type{1};
  constexpr auto sourceDof = size_type{2};
  constexpr auto value = value_type{.7};
  constexpr auto factor = value_type{2.};

  const auto result = boundaryConditionsToTransform(
      {
          {{targetVertex, targetDof},
           {{factor, {sourceVertex, sourceDof}}},
           0.},
      },
      this->mesh);

  const auto image = this->mesh.toCanonicalOrder(apply(result, [&]() {
    auto input = createTransformInput(result.matrix);
    input.unwrap().fill(value);
    return input;
  }()));

  EXPECT_THAT(
      image.unwrap(),
      ScalarEqIfLocal(sourceVertex * this->dofPerVertex + sourceDof, value));
  EXPECT_THAT(image.unwrap(),
              ScalarEqIfLocal(targetVertex * this->dofPerVertex + targetDof,
                              factor * value));
}

TYPED_TEST(boundaryConditionsToTransform_Test,
           dof_can_be_set_to_affine_combination_of_other_dofs) {
  using size_type = typename TestFixture::size_type;
  using value_type = typename TestFixture::value_type;

  constexpr auto targetVertex = size_type{0};
  constexpr auto targetDof = size_type{1};
  constexpr auto factor_1 = value_type{2.};
  constexpr auto sourceVertex_1 = size_type{1};
  constexpr auto sourceDof_1 = size_type{1};
  constexpr auto factor_2 = value_type{3.};
  constexpr auto sourceVertex_2 = size_type{1};
  constexpr auto sourceDof_2 = size_type{2};
  constexpr auto offset = value_type{.3};
  constexpr auto value = value_type{.7};

  const auto result = boundaryConditionsToTransform(
      {
          {{targetVertex, targetDof},
           {
               {factor_1, {sourceVertex_1, sourceDof_1}},
               {factor_2, {sourceVertex_2, sourceDof_2}},
           },
           offset},
      },
      this->mesh);

  const auto image = this->mesh.toCanonicalOrder(apply(result, [&]() {
    auto input = createTransformInput(result.matrix);
    input.unwrap().fill(value);
    return input;
  }()));

  EXPECT_THAT(image.unwrap(),
              ScalarEqIfLocal(sourceVertex_1 * this->dofPerVertex + sourceDof_1,
                              value));
  EXPECT_THAT(image.unwrap(),
              ScalarEqIfLocal(sourceVertex_2 * this->dofPerVertex + sourceDof_2,
                              value));
  EXPECT_THAT(image.unwrap(),
              ScalarEqIfLocal(targetVertex * this->dofPerVertex + targetDof,
                              factor_1 * value + factor_2 * value + offset));
}

TYPED_TEST(boundaryConditionsToTransform_Test,
           throws_if_source_is_also_target) {
  using size_type = typename TestFixture::size_type;
  using value_type = typename TestFixture::value_type;

  constexpr auto targetVertex = size_type{0};
  constexpr auto targetDof = size_type{1};
  constexpr auto sourceVertex = size_type{1};
  constexpr auto sourceDof = size_type{2};
  constexpr auto value = value_type{.7};

  const auto call = [&]() {
    return boundaryConditionsToTransform(
        {
            {{targetVertex, targetDof},
             {{1., {sourceVertex, sourceDof}}},
             value},
            {{sourceVertex, sourceDof}, {}, value},
        },
        this->mesh);
  };

  EXPECT_THROW(call(), InconsistentBoundaryConditionsException);
}

TYPED_TEST(boundaryConditionsToTransform_Test,
           throws_if_target_vertex_is_invalid) {
  using size_type = typename TestFixture::size_type;
  using value_type = typename TestFixture::value_type;

  constexpr auto targetVertex = size_type{7};
  constexpr auto targetDof = size_type{1};
  constexpr auto value = value_type{.7};

  const auto call = [&]() {
    return boundaryConditionsToTransform(
        {
            {{targetVertex, targetDof}, {}, value},
        },
        this->mesh);
  };

  EXPECT_THROW(call(), InvalidVertexException);
}

TYPED_TEST(boundaryConditionsToTransform_Test,
           throws_if_source_vertex_is_invalid) {
  using size_type = typename TestFixture::size_type;

  constexpr auto targetVertex = size_type{0};
  constexpr auto targetDof = size_type{1};
  constexpr auto sourceVertex = size_type{7};
  constexpr auto sourceDof = size_type{2};

  const auto call = [&]() {
    return boundaryConditionsToTransform(
        {
            {{targetVertex, targetDof}, {{1., {sourceVertex, sourceDof}}}, 0.},
        },
        this->mesh);
  };

  EXPECT_THROW(call(), InvalidVertexException);
}

} // namespace
} // namespace solve
} // namespace ae108