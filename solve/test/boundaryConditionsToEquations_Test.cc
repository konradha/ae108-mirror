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
#include "ae108/solve/InvalidVertexException.h"
#include "ae108/solve/boundaryConditionsToEquations.h"
#include <array>
#include <gmock/gmock.h>

using ae108::cppptest::ScalarEqIfLocal;
using testing::DoubleEq;
using testing::Eq;
using testing::Ge;
using testing::Pair;
using testing::Test;

namespace ae108 {
namespace solve {
namespace {

template <class Policy> struct boundaryConditionsToEquations_Test : Test {
  using mesh_type = cpppetsc::Mesh<Policy>;
  using size_type = typename mesh_type::size_type;
  using value_type = typename mesh_type::value_type;
  using vector_type = typename mesh_type::vector_type;

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

TYPED_TEST_CASE(boundaryConditionsToEquations_Test, Policies);

TYPED_TEST(boundaryConditionsToEquations_Test,
           returns_size_zero_matrix_if_no_conditions) {
  const auto result = boundaryConditionsToEquations({}, this->mesh);

  EXPECT_THAT(result.matrix.size(),
              Pair(0, this->numberOfVertices * this->dofPerVertex));
}

TYPED_TEST(boundaryConditionsToEquations_Test,
           returns_size_zero_vector_if_no_conditions) {
  const auto result = boundaryConditionsToEquations({}, this->mesh);

  EXPECT_THAT(result.vector.unwrap().size(), Eq(0));
}

TYPED_TEST(boundaryConditionsToEquations_Test,
           returns_one_row_per_rank_matrix_if_one_condition) {
  using size_type = typename TestFixture::size_type;
  using value_type = typename TestFixture::value_type;

  constexpr auto vertex = size_type{1};
  constexpr auto dof = size_type{2};
  constexpr auto value = value_type{.7};

  const auto result =
      boundaryConditionsToEquations({{{vertex, dof}, {}, value}}, this->mesh);

  EXPECT_THAT(result.matrix.localSize().first, Eq(1));
}

TYPED_TEST(boundaryConditionsToEquations_Test,
           returns_one_row_per_rank_vector_if_one_condition) {
  using size_type = typename TestFixture::size_type;
  using value_type = typename TestFixture::value_type;

  constexpr auto vertex = size_type{1};
  constexpr auto dof = size_type{2};
  constexpr auto value = value_type{.7};

  const auto result =
      boundaryConditionsToEquations({{{vertex, dof}, {}, value}}, this->mesh);

  EXPECT_THAT(result.vector.unwrap().localSize(), Eq(1));
}

TYPED_TEST(boundaryConditionsToEquations_Test,
           throws_invalid_value_exception_if_invalid_target) {
  using size_type = typename TestFixture::size_type;
  using value_type = typename TestFixture::value_type;

  constexpr auto vertex = size_type{1234};
  constexpr auto dof = size_type{2};
  constexpr auto value = value_type{.7};

  EXPECT_THROW(
      boundaryConditionsToEquations({{{vertex, dof}, {}, value}}, this->mesh),
      InvalidVertexException);
}

TYPED_TEST(boundaryConditionsToEquations_Test,
           application_of_transform_to_zero_input_yields_offset) {
  using size_type = typename TestFixture::size_type;
  using value_type = typename TestFixture::value_type;

  constexpr auto vertex = size_type{1};
  constexpr auto dof = size_type{2};
  constexpr auto value = value_type{.7};

  const auto result =
      boundaryConditionsToEquations({{{vertex, dof}, {}, value}}, this->mesh);

  const auto image = apply(result, createTransformInput(result.matrix));
  for (auto row = size_type{0}; row < image.unwrap().size(); ++row) {
    EXPECT_THAT(image.unwrap(), ScalarEqIfLocal(row, -1. * value));
  }
}

TYPED_TEST(boundaryConditionsToEquations_Test,
           application_of_transform_yields_offset_for_multiple_conditions) {
  using size_type = typename TestFixture::size_type;
  using value_type = typename TestFixture::value_type;
  using vector_type = typename TestFixture::vector_type;

  constexpr std::array<size_type, 2> vertices = {{0, 1}};
  constexpr std::array<size_type, 2> dofs = {{1, 2}};
  constexpr std::array<value_type, 2> values = {{.7, .8}};
  constexpr auto factor = .4;

  auto displacements = vector_type::fromGlobalMesh(this->mesh);
  auto localDisplacements = vector_type::fromLocalMesh(this->mesh);
  auto boundaryConditions =
      std::vector<cpppetsc::GeneralizedMeshBoundaryCondition<
          typename TestFixture::mesh_type>>{};

  for (auto &&vertex : this->mesh.localVertices()) {
    switch (vertex.index()) {
    case vertices[0]:
      boundaryConditions.push_back({{vertices[0], dofs[0]}, {}, values[0]});
      boundaryConditions.push_back({{vertices[0], dofs[1]},
                                    {{factor, {vertices[1], dofs[1]}}},
                                    values[1]});
      vertex.setVertexData({0., values[0], .3 * factor + values[1]},
                           &localDisplacements);
      break;
    case vertices[1]:
      vertex.setVertexData({.1, .2, .3}, &localDisplacements);
      break;
    }
  }
  this->mesh.copyToGlobalVector(localDisplacements, &displacements);

  const auto result =
      boundaryConditionsToEquations(boundaryConditions, this->mesh);
  const auto image = apply(result, displacements);

  for (auto row = size_type{0}; row < image.unwrap().size(); ++row) {
    EXPECT_THAT(image.unwrap(), ScalarEqIfLocal(row, 0.));
  }
}

} // namespace
} // namespace solve
} // namespace ae108