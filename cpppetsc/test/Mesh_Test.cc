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

#include "ae108/cpppetsc/Matrix.h"
#include "ae108/cpppetsc/Mesh.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cppptest/Matchers.h"
#include <gmock/gmock.h>
#include <numeric>
#include <vector>

using ae108::cppptest::AlmostEqIfLocal;
using testing::DoubleEq;
using testing::Each;
using testing::ElementsAre;
using testing::ElementsAreArray;
using testing::Eq;
using testing::Ge;
using testing::Le;
using testing::Not;
using testing::Pair;
using testing::Pointwise;
using testing::SizeIs;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace cpppetsc {
namespace {

/**
 * @brief Create a vector of the given length with contents:
 * increment, 2 * increment, ...
 */
template <class ValueType>
std::vector<ValueType>
createData(typename std::vector<ValueType>::size_type size,
           ValueType increment = ValueType{1.}) {
  std::vector<ValueType> result(size, increment);
  std::partial_sum(result.begin(), result.end(), result.begin());
  return result;
}

/**
 * @brief Create a distributed vector that contains the vertex/element index
 * for every vertex/element dof.
 */
template <class Policy>
distributed<typename Mesh<Policy>::vector_type>
createIndexVector(const Mesh<Policy> &mesh) {
  using vector_type = typename Mesh<Policy>::vector_type;
  using value_type = typename Mesh<Policy>::value_type;

  auto vector = vector_type::fromLocalMesh(mesh);
  for (const auto &vertex : mesh.localVertices()) {
    const auto values = std::vector<value_type>(
        vertex.numberOfDofs(), static_cast<value_type>(vertex.index()));
    vertex.setVertexData(values, &vector);
  }
  for (const auto &element : mesh.localElements()) {
    const auto range = element.localDofLineRange();
    const auto replacer = vector.unwrap().replace();
    for (auto i = range.first; i < range.second; ++i) {
      replacer(i) = static_cast<value_type>(element.index());
    }
  }

  auto result = vector_type::fromGlobalMesh(mesh);
  mesh.copyToGlobalVector(vector, &result);
  return result;
}

/**
 * @brief Returns a vector of the local values in parameter.
 */
template <class Policy>
std::vector<typename Vector<Policy>::value_type>
localValues(const Vector<Policy> &input) {
  return std::vector<typename Vector<Policy>::value_type>(input.localBegin(),
                                                          input.localEnd());
}

template <typename T> struct Mesh_Test : Test {
  using policy_type = typename T::type;
  using mesh_type = Mesh<policy_type>;

  using vector_type = typename mesh_type::vector_type;
  using matrix_type = typename mesh_type::matrix_type;
  using size_type = typename mesh_type::size_type;
  using value_type = typename mesh_type::value_type;

  static constexpr size_type dimension = 1;
  static constexpr size_type totalNumberOfElements = 3;
  static constexpr size_type totalNumberOfVertices = 3;
  static constexpr size_type verticesPerElement = 2;
  static constexpr size_type dofPerVertex = T::DofPerVertex;
  static constexpr size_type dofPerElement = T::DofPerElement;

  using Connectivity = std::vector<std::array<size_type, verticesPerElement>>;
  const Connectivity connectivity = {{0, 2}, {1, 2}, {0, 2}};

  const mesh_type mesh = mesh_type::fromConnectivity(
      dimension, connectivity, totalNumberOfVertices, dofPerVertex,
      dofPerElement);
};

template <class Policy, int VertexDof, int ElementDof> struct TestCase {
  using type = Policy;
  static constexpr auto DofPerVertex = VertexDof;
  static constexpr auto DofPerElement = ElementDof;
};

using TestCases = Types<TestCase<SequentialComputePolicy, 2, 0>,
                        TestCase<SequentialComputePolicy, 0, 3>,
                        TestCase<SequentialComputePolicy, 2, 3>,
                        TestCase<ParallelComputePolicy, 2, 0>,
                        TestCase<ParallelComputePolicy, 0, 3>,
                        TestCase<ParallelComputePolicy, 2, 3>>;
TYPED_TEST_CASE(Mesh_Test, TestCases);

TYPED_TEST(Mesh_Test, from_connectivity_generates_correct_number_of_entities) {
  EXPECT_THAT(this->mesh.totalNumberOfElements(),
              Eq(this->totalNumberOfElements));
  EXPECT_THAT(this->mesh.totalNumberOfVertices(),
              Eq(this->totalNumberOfVertices));
}

TYPED_TEST(Mesh_Test, cloned_mesh_has_correct_number_of_entities) {
  using size_type = typename TestFixture::size_type;

  const auto mesh = this->mesh.cloneWithDofs(size_type{0}, size_type{0});
  EXPECT_THAT(mesh.totalNumberOfElements(), Eq(this->totalNumberOfElements));
  EXPECT_THAT(mesh.totalNumberOfVertices(), Eq(this->totalNumberOfVertices));
}

TYPED_TEST(Mesh_Test, cloned_mesh_has_correct_connectivity) {
  using size_type = typename TestFixture::size_type;

  const auto clonedMesh = this->mesh.cloneWithDofs(size_type{0}, size_type{0});
  const auto clonedRange = clonedMesh.localElements();
  auto iterator = clonedRange.begin();

  for (const auto &element : this->mesh.localElements()) {
    ASSERT_THAT(iterator, Not(Eq(clonedRange.end())));

    EXPECT_THAT(iterator->index(), Eq(element.index()));
    EXPECT_THAT(iterator->vertexIndices(), Eq(element.vertexIndices()));

    ++iterator;
  }
}

TYPED_TEST(Mesh_Test, view_for_cloned_mesh_returns_same_vertex_index) {
  using size_type = typename TestFixture::size_type;

  const auto clonedMesh = this->mesh.cloneWithDofs(size_type{0}, size_type{0});

  for (const auto &vertex : this->mesh.localVertices()) {
    EXPECT_THAT(vertex.forClonedMesh(&clonedMesh).index(), Eq(vertex.index()));
  }
}

TYPED_TEST(Mesh_Test, correct_number_of_local_elements_are_iterated) {
  using size_type = typename TestFixture::size_type;

  const auto range = this->mesh.localElements();
  auto elements = size_type{std::distance(range.begin(), range.end())};

  ASSERT_THAT(MPI_Allreduce(MPI_IN_PLACE, &elements, 1, MPI_INT, MPI_SUM,
                            TestFixture::policy_type::communicator()),
              Eq(0));
  EXPECT_THAT(elements, Eq(this->totalNumberOfElements));
}

TYPED_TEST(Mesh_Test, plausible_number_of_local_vertices_are_iterated) {
  using size_type = typename TestFixture::size_type;

  const auto range = this->mesh.localVertices();
  auto vertices = size_type{std::distance(range.begin(), range.end())};

  ASSERT_THAT(MPI_Allreduce(MPI_IN_PLACE, &vertices, 1, MPI_INT, MPI_SUM,
                            TestFixture::policy_type::communicator()),
              Eq(0));
  EXPECT_THAT(vertices, Ge(this->totalNumberOfVertices));
}

TYPED_TEST(Mesh_Test, vertices_have_correct_number_of_dofs) {
  for (const auto &vertex : this->mesh.localVertices()) {
    EXPECT_THAT(vertex.numberOfDofs(), Eq(this->dofPerVertex));
  }
}

TYPED_TEST(Mesh_Test, clone_has_correct_vertex_dofs) {
  using size_type = typename TestFixture::size_type;

  const auto dofPerVertex = size_type{77};
  const auto dofPerElement = size_type{777};
  const auto mesh = this->mesh.cloneWithDofs(dofPerVertex, dofPerElement);
  for (const auto &vertex : mesh.localVertices()) {
    EXPECT_THAT(vertex.numberOfDofs(), Eq(dofPerVertex));
  }
}

TYPED_TEST(Mesh_Test, elements_have_correct_number_of_dofs) {
  for (const auto &element : this->mesh.localElements()) {
    EXPECT_THAT(element.numberOfDofs(), Eq(this->dofPerElement));
  }
}

TYPED_TEST(Mesh_Test, clone_has_correct_element_dofs) {
  using size_type = typename TestFixture::size_type;

  const auto dofPerVertex = size_type{77};
  const auto dofPerElement = size_type{777};
  const auto mesh = this->mesh.cloneWithDofs(dofPerVertex, dofPerElement);
  for (const auto &element : mesh.localElements()) {
    EXPECT_THAT(element.numberOfDofs(), Eq(dofPerElement));
  }
}

TYPED_TEST(Mesh_Test, local_number_of_vertices_is_plausible) {
  auto result = this->mesh.localNumberOfVertices();
  ASSERT_THAT(MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_INT, MPI_SUM,
                            TestFixture::policy_type::communicator()),
              Eq(0));

  EXPECT_THAT(result, Ge(this->totalNumberOfVertices));
}

TYPED_TEST(Mesh_Test, local_number_of_elements_is_correct) {
  auto result = this->mesh.localNumberOfElements();
  ASSERT_THAT(MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_INT, MPI_SUM,
                            TestFixture::policy_type::communicator()),
              Eq(0));

  EXPECT_THAT(result, Eq(this->totalNumberOfElements));
}

TYPED_TEST(Mesh_Test, created_matrix_has_correct_size) {
  using matrix_type = typename TestFixture::matrix_type;
  auto matrix = matrix_type::fromMesh(this->mesh);

  const auto reference = this->totalNumberOfVertices * this->dofPerVertex +
                         this->totalNumberOfElements * this->dofPerElement;
  EXPECT_THAT(matrix.size(), Pair(Eq(reference), Eq(reference)));
}

TYPED_TEST(Mesh_Test, created_global_vector_has_correct_size) {
  using vector_type = typename TestFixture::vector_type;
  const auto vector = vector_type::fromGlobalMesh(this->mesh);

  const auto reference = this->totalNumberOfVertices * this->dofPerVertex +
                         this->totalNumberOfElements * this->dofPerElement;
  EXPECT_THAT(vector.unwrap(), SizeIs(reference));
}

TYPED_TEST(Mesh_Test, global_vectors_are_initialized_to_zero) {
  using vector_type = typename TestFixture::vector_type;
  const auto vector = vector_type::fromGlobalMesh(this->mesh);

  EXPECT_THAT(localValues(vector.unwrap()), Each(DoubleEq(0.)));
}

TYPED_TEST(Mesh_Test, value_can_be_read_from_global_vector) {
  using vector_type = typename TestFixture::vector_type;
  using value_type = typename TestFixture::value_type;

  auto vector = vector_type::fromGlobalMesh(this->mesh);
  const auto value = value_type{.7};
  vector.unwrap().fill(value);

  const auto fullVector = vector_type::fromDistributedInCanonicalOrder(
      vector, this->mesh, this->dofPerVertex, this->dofPerElement);

  const auto reference = this->totalNumberOfVertices * this->dofPerVertex +
                         this->totalNumberOfElements * this->dofPerElement;
  EXPECT_THAT(localValues(fullVector.unwrap()), SizeIs(reference));
  EXPECT_THAT(localValues(fullVector.unwrap()), Each(DoubleEq(value)));
}

TYPED_TEST(Mesh_Test, created_local_vector_has_plausible_size) {
  using vector_type = typename TestFixture::vector_type;
  using size_type = typename TestFixture::size_type;
  auto size = vector_type::fromLocalMesh(this->mesh).unwrap().size();

  auto sum = size_type{0};

  ASSERT_THAT(MPI_Allreduce(&size, &sum, 1, MPI_INT, MPI_SUM,
                            TestFixture::policy_type::communicator()),
              Eq(0));

  const auto reference = this->totalNumberOfVertices * this->dofPerVertex +
                         this->totalNumberOfElements * this->dofPerElement;
  EXPECT_THAT(sum, Ge(reference));
}

TYPED_TEST(Mesh_Test, local_vectors_are_initialized_to_zero) {
  using vector_type = typename TestFixture::vector_type;
  const auto vector = vector_type::fromLocalMesh(this->mesh);

  EXPECT_THAT(localValues(vector.unwrap()), Each(DoubleEq(0.)));
}

TYPED_TEST(Mesh_Test, value_can_be_read_from_local_vector) {
  using vector_type = typename TestFixture::vector_type;
  using value_type = typename TestFixture::value_type;

  auto vector = vector_type::fromLocalMesh(this->mesh);
  const auto value = value_type{.7};
  vector.unwrap().fill(value);

  const auto reference = this->totalNumberOfVertices * this->dofPerVertex +
                         this->totalNumberOfElements * this->dofPerElement;
  EXPECT_THAT(localValues(vector.unwrap()), SizeIs(Le(reference)));
  EXPECT_THAT(localValues(vector.unwrap()), Each(DoubleEq(value)));
}

TYPED_TEST(Mesh_Test, local_entities_correspond_to_local_vector_size) {
  using vector_type = typename TestFixture::vector_type;

  const auto vector = vector_type::fromLocalMesh(this->mesh);

  const auto reference =
      this->mesh.localNumberOfVertices() * this->dofPerVertex +
      this->mesh.localNumberOfElements() * this->dofPerElement;
  EXPECT_THAT(vector.unwrap().size(), Eq(reference));
}

TYPED_TEST(Mesh_Test, copying_to_local_vector_works) {
  using vector_type = typename TestFixture::vector_type;
  using value_type = typename TestFixture::value_type;

  auto globalVector = vector_type::fromGlobalMesh(this->mesh);
  const auto value = value_type{.7};
  globalVector.unwrap().fill(value);

  auto localVector = vector_type::fromLocalMesh(this->mesh);
  this->mesh.copyToLocalVector(globalVector, &localVector);

  EXPECT_THAT(localValues(localVector.unwrap()), Each(DoubleEq(value)));
}

TYPED_TEST(Mesh_Test, adding_to_global_vector_works) {
  using vector_type = typename TestFixture::vector_type;
  using size_type = typename TestFixture::size_type;
  using value_type = typename TestFixture::value_type;

  auto localVector = vector_type::fromLocalMesh(this->mesh);
  const auto value = value_type{.7};
  localVector.unwrap().fill(value);

  auto globalVector = vector_type::fromGlobalMesh(this->mesh);
  const auto value_2 = value_type{.3};
  globalVector.unwrap().fill(value_2);

  this->mesh.addToGlobalVector(localVector, &globalVector);

  const auto fullVector = vector_type::fromDistributedInCanonicalOrder(
      globalVector, this->mesh, this->dofPerVertex, this->dofPerElement);

  ASSERT_THAT(fullVector.unwrap(),
              SizeIs(this->totalNumberOfVertices * this->dofPerVertex +
                     this->totalNumberOfElements * this->dofPerElement));

  auto index = size_type{0};
  for (size_type i = 0; i < this->dofPerVertex; ++i) {
    EXPECT_THAT(fullVector(index++), Ge(value + value_2));
  }
  for (size_type i = 0; i < this->dofPerVertex; ++i) {
    EXPECT_THAT(fullVector(index++), DoubleEq(value + value_2));
  }
  for (size_type i = 0; i < this->dofPerVertex; ++i) {
    EXPECT_THAT(fullVector(index++), Ge(value + value_2));
  }
  for (size_type i = 0; i < this->dofPerElement; ++i) {
    EXPECT_THAT(fullVector(index++), DoubleEq(value + value_2));
  }
}

TYPED_TEST(Mesh_Test, copying_to_global_vector_works) {
  using vector_type = typename TestFixture::vector_type;
  using value_type = typename TestFixture::value_type;

  auto localVector = vector_type::fromLocalMesh(this->mesh);
  const auto value = value_type{.7};
  localVector.unwrap().fill(value);

  auto globalVector = vector_type::fromGlobalMesh(this->mesh);
  const auto value_2 = value_type{.3};
  globalVector.unwrap().fill(value_2);

  this->mesh.copyToGlobalVector(localVector, &globalVector);

  const auto fullVector = vector_type::fromDistributedInCanonicalOrder(
      globalVector, this->mesh, this->dofPerVertex, this->dofPerElement);

  ASSERT_THAT(fullVector.unwrap(),
              SizeIs(this->totalNumberOfVertices * this->dofPerVertex +
                     this->totalNumberOfElements * this->dofPerElement));
  EXPECT_THAT(localValues(fullVector.unwrap()), Each(DoubleEq(value)));
}

TYPED_TEST(Mesh_Test, begin_element_iterator_starts_at_index_0) {
  using size_type = typename TestFixture::size_type;
  auto iterator = this->mesh.localElements().begin();

  auto result = std::numeric_limits<size_type>::max();

  if (iterator != this->mesh.localElements().end()) {
    result = iterator->index();
  }
  ASSERT_THAT(MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_INT, MPI_MIN,
                            TestFixture::policy_type::communicator()),
              Eq(0));

  EXPECT_THAT(result, Eq(0));
}

TYPED_TEST(Mesh_Test, end_element_iterator_is_past_final_index) {
  using size_type = typename TestFixture::size_type;
  auto iterator = this->mesh.localElements().end();

  auto result = std::numeric_limits<size_type>::min();

  if (iterator != this->mesh.localElements().begin()) {
    result = (--iterator)->index();
  }
  ASSERT_THAT(MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_INT, MPI_MAX,
                            TestFixture::policy_type::communicator()),
              Eq(0));

  EXPECT_THAT(result, Eq(this->totalNumberOfElements - 1));
}

TYPED_TEST(Mesh_Test, begin_vertex_iterator_starts_at_index_0) {
  using size_type = typename TestFixture::size_type;
  auto iterator = this->mesh.localVertices().begin();

  auto result = std::numeric_limits<size_type>::max();

  if (iterator != this->mesh.localVertices().end()) {
    result = iterator->index();
  }
  ASSERT_THAT(MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_INT, MPI_MIN,
                            TestFixture::policy_type::communicator()),
              Eq(0));

  EXPECT_THAT(result, Eq(0));
}

TYPED_TEST(Mesh_Test, end_vertex_iterator_is_past_final_index) {
  using size_type = typename TestFixture::size_type;
  auto iterator = this->mesh.localVertices().end();

  auto result = std::numeric_limits<size_type>::min();

  if (iterator != this->mesh.localVertices().begin()) {
    result = (--iterator)->index();
  }
  ASSERT_THAT(MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_INT, MPI_MAX,
                            TestFixture::policy_type::communicator()),
              Eq(0));

  EXPECT_THAT(result, Eq(this->totalNumberOfVertices - 1));
}

TYPED_TEST(Mesh_Test, vertex_dofs_add_up_to_ge_problem_size) {
  using size_type = typename TestFixture::size_type;

  auto result = size_type{0};

  for (const auto &vertex : this->mesh.localVertices()) {
    result += vertex.numberOfDofs();
  }

  ASSERT_THAT(MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_INT, MPI_SUM,
                            TestFixture::policy_type::communicator()),
              Eq(0));

  const auto problemSize = this->dofPerVertex * this->totalNumberOfVertices;
  EXPECT_THAT(result, Ge(problemSize));
}

TYPED_TEST(Mesh_Test, element_dofs_add_up_to_problem_size) {
  using size_type = typename TestFixture::size_type;

  auto result = size_type{0};

  for (const auto &element : this->mesh.localElements()) {
    result += element.numberOfDofs();
  }

  ASSERT_THAT(MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_INT, MPI_SUM,
                            TestFixture::policy_type::communicator()),
              Eq(0));

  const auto problemSize = this->dofPerElement * this->totalNumberOfElements;
  EXPECT_THAT(result, Eq(problemSize));
}

TYPED_TEST(Mesh_Test, local_line_ranges_have_correct_size) {
  for (const auto &vertex : this->mesh.localVertices()) {
    const auto range = vertex.localDofLineRange();
    EXPECT_THAT(range.second - range.first, Eq(vertex.numberOfDofs()));
  }
  for (const auto &element : this->mesh.localElements()) {
    const auto range = element.localDofLineRange();
    EXPECT_THAT(range.second - range.first, Eq(element.numberOfDofs()));
  }
}

TYPED_TEST(Mesh_Test, global_line_ranges_have_correct_size) {
  for (const auto &vertex : this->mesh.localVertices()) {
    const auto range = vertex.globalDofLineRange();
    EXPECT_THAT(std::abs(range.second - range.first),
                Eq(vertex.numberOfDofs()));
  }
  for (const auto &element : this->mesh.localElements()) {
    const auto range = element.globalDofLineRange();
    EXPECT_THAT(std::abs(range.second - range.first),
                Eq(element.numberOfDofs()));
  }
}

TYPED_TEST(Mesh_Test, local_line_range_accesses_correct_region) {
  const auto globalVector = createIndexVector(this->mesh);

  using vector_type = typename TestFixture::vector_type;
  auto localVector = vector_type::fromLocalMesh(this->mesh);

  this->mesh.copyToLocalVector(globalVector, &localVector);

  for (const auto &vertex : this->mesh.localVertices()) {
    const auto range = vertex.localDofLineRange();
    for (auto i = range.first; i < range.second; ++i) {
      EXPECT_THAT(localVector(i), DoubleEq(vertex.index()));
    }
  }
  for (const auto &element : this->mesh.localElements()) {
    const auto range = element.localDofLineRange();
    for (auto i = range.first; i < range.second; ++i) {
      EXPECT_THAT(localVector(i), DoubleEq(element.index()));
    }
  }
}

TYPED_TEST(Mesh_Test, global_line_range_accesses_correct_region) {
  const auto globalVector = createIndexVector(this->mesh);

  for (const auto &vertex : this->mesh.localVertices()) {
    const auto range = vertex.globalDofLineRange();
    for (auto i = range.first; i < range.second; ++i) {
      EXPECT_THAT(globalVector(i), DoubleEq(vertex.index()));
    }
  }
  for (const auto &element : this->mesh.localElements()) {
    const auto range = element.globalDofLineRange();
    for (auto i = range.first; i < range.second; ++i) {
      EXPECT_THAT(globalVector(i), DoubleEq(element.index()));
    }
  }
}

TYPED_TEST(Mesh_Test, sorting_in_canonical_order_works) {
  using vector_type = typename TestFixture::vector_type;
  using size_type = typename TestFixture::size_type;

  const auto globalVector = createIndexVector(this->mesh);
  const auto fullVector = vector_type::fromDistributedInCanonicalOrder(
      globalVector, this->mesh, this->dofPerVertex, this->dofPerElement);

  ASSERT_THAT(fullVector.unwrap(),
              SizeIs(this->totalNumberOfVertices * this->dofPerVertex +
                     this->totalNumberOfElements * this->dofPerElement));
  for (size_type i = 0; i < this->totalNumberOfVertices * this->dofPerVertex;
       ++i) {
    EXPECT_THAT(fullVector(i),
                DoubleEq(static_cast<double>(i / this->dofPerVertex)));
  }
  for (size_type i = 0; i < this->totalNumberOfElements * this->dofPerElement;
       ++i) {
    EXPECT_THAT(
        fullVector(this->totalNumberOfVertices * this->dofPerVertex + i),
        DoubleEq(static_cast<double>(i / this->dofPerElement)));
  }
}

TYPED_TEST(Mesh_Test, vertex_indices_are_correct_for_index) {
  for (const auto &element : this->mesh.localElements()) {
    EXPECT_THAT(element.vertexIndices(),
                ElementsAreArray(this->connectivity.at(element.index())));
  }
}

TYPED_TEST(Mesh_Test, number_of_vertices_is_correct_per_element) {
  for (const auto &element : this->mesh.localElements()) {
    EXPECT_THAT(element.numberOfVertices(), Eq(this->verticesPerElement));
    EXPECT_THAT(element.vertexIndices().size(), Eq(element.numberOfVertices()));
  }
}

TYPED_TEST(Mesh_Test, add_transferring_element_data_from_and_to_vector_works) {
  using vector_type = typename TestFixture::vector_type;
  using value_type = typename TestFixture::value_type;

  for (const auto &element : this->mesh.localElements()) {
    auto localVector = vector_type::fromLocalMesh(this->mesh);

    const auto input = createData<value_type>(
        this->verticesPerElement * this->dofPerVertex + this->dofPerElement);
    element.addElementData(input, &localVector);
    element.addElementData(input, &localVector);

    std::vector<value_type> output;
    element.copyElementData(localVector, &output);

    const auto reference = createData<value_type>(
        this->verticesPerElement * this->dofPerVertex + this->dofPerElement,
        value_type{2.});
    EXPECT_THAT(output, Pointwise(DoubleEq(), reference));
  }
}

TYPED_TEST(Mesh_Test, set_transferring_element_data_from_and_to_vector_works) {
  using vector_type = typename TestFixture::vector_type;
  using value_type = typename TestFixture::value_type;

  for (const auto &element : this->mesh.localElements()) {
    auto localVector = vector_type::fromLocalMesh(this->mesh);

    const auto input = createData<value_type>(
        this->verticesPerElement * this->dofPerVertex + this->dofPerElement);
    element.setElementData(input, &localVector);

    std::vector<value_type> output;
    element.copyElementData(localVector, &output);

    EXPECT_THAT(output, Eq(input));
  }
}

TYPED_TEST(Mesh_Test, add_transferring_vertex_data_from_and_to_vector_works) {
  using vector_type = typename TestFixture::vector_type;
  using value_type = typename TestFixture::value_type;

  for (const auto &vertex : this->mesh.localVertices()) {
    auto localVector = vector_type::fromLocalMesh(this->mesh);

    const auto input = createData<value_type>(this->dofPerVertex);
    vertex.addVertexData(input, &localVector);
    vertex.addVertexData(input, &localVector);

    std::vector<value_type> output;
    vertex.copyVertexData(localVector, &output);

    const auto reference =
        createData<value_type>(this->dofPerVertex, value_type{2.});
    EXPECT_THAT(output, Pointwise(DoubleEq(), reference));
  }
}

TYPED_TEST(Mesh_Test, set_transferring_vertex_data_from_and_to_vector_works) {
  using vector_type = typename TestFixture::vector_type;
  using value_type = typename TestFixture::value_type;

  for (const auto &vertex : this->mesh.localVertices()) {
    auto localVector = vector_type::fromLocalMesh(this->mesh);

    const auto input = createData<value_type>(this->dofPerVertex);
    vertex.setVertexData(input, &localVector);

    std::vector<value_type> output;
    vertex.copyVertexData(localVector, &output);

    EXPECT_THAT(output, Eq(input));
  }
}

TYPED_TEST(Mesh_Test, adding_element_data_targets_local_dof_range) {
  using vector_type = typename TestFixture::vector_type;
  using value_type = typename TestFixture::value_type;

  for (const auto &element : this->mesh.localElements()) {
    auto localVector = vector_type::fromLocalMesh(this->mesh);
    const auto input = createData<value_type>(
        this->verticesPerElement * this->dofPerVertex + this->dofPerElement);
    element.addElementData(input, &localVector);

    const auto range = element.localDofLineRange();
    for (auto index = range.first; index < range.second; ++index) {
      EXPECT_THAT(localVector(index), DoubleEq(input.at(index - range.first)));
    }
  }
}

TYPED_TEST(Mesh_Test, adding_vertex_data_target_local_dof_range) {
  using vector_type = typename TestFixture::vector_type;
  using value_type = typename TestFixture::value_type;

  for (const auto &vertex : this->mesh.localVertices()) {
    auto localVector = vector_type::fromLocalMesh(this->mesh);
    const auto input = createData<value_type>(this->dofPerVertex);
    vertex.addVertexData(input, &localVector);

    const auto range = vertex.localDofLineRange();
    for (auto index = range.first; index < range.second; ++index) {
      EXPECT_THAT(localVector(index), DoubleEq(input.at(index - range.first)));
    }
  }
}

TYPED_TEST(Mesh_Test, adding_element_0_to_matrix_works) {
  using value_type = typename TestFixture::value_type;
  using size_type = typename TestFixture::size_type;
  using matrix_type = typename TestFixture::matrix_type;

  const auto onlyVertexDofs = bool{this->dofPerElement == size_type{0}};
  if (!onlyVertexDofs) {
    return;
  }

  const auto elementIndex = size_type{0};
  auto result = matrix_type::fromMesh(this->mesh);

  for (const auto &element : this->mesh.localElements()) {
    if (element.index() != elementIndex) {
      continue;
    }

    ASSERT_THAT(element.vertexIndices(),
                ElementsAre(this->connectivity[elementIndex][0],
                            this->connectivity[elementIndex][1]));
    const auto input = createData<value_type>(
        (this->verticesPerElement * this->dofPerVertex + this->dofPerElement) *
        (this->verticesPerElement * this->dofPerVertex + this->dofPerElement));
    element.addElementMatrix(input, &result);
  }

  result.finalize();

  const auto vertexIndex_1 = this->connectivity[elementIndex][0];
  const auto vertexIndex_2 = this->connectivity[elementIndex][1];

  for (const auto &vertex : this->mesh.localVertices()) {
    if (vertex.index() != vertexIndex_1 && vertex.index() != vertexIndex_2) {
      continue;
    }

    const auto range = vertex.globalDofLineRange();
    auto reference = value_type{vertex.index() == vertexIndex_1 ? 1. : 11.};

    for (auto i = range.first; i < range.second; ++i) {
      for (auto j = range.first; j < range.second; ++j) {
        EXPECT_THAT(result, AlmostEqIfLocal(i, j, reference));
        reference++;
      }
      reference += 2.;
    }
  }
}

TYPED_TEST(Mesh_Test, adding_vertex_0_to_matrix_works) {
  using value_type = typename TestFixture::value_type;
  using matrix_type = typename TestFixture::matrix_type;

  const auto vertexIndex = 0;
  auto result = matrix_type::fromMesh(this->mesh);

  for (const auto &vertex : this->mesh.localVertices()) {
    const bool isGhost(vertex.globalDofLineRange().first < 0);
    if (vertex.index() != vertexIndex || isGhost) {
      continue;
    }

    const auto input =
        createData<value_type>(this->verticesPerElement * this->dofPerVertex);
    vertex.addVertexMatrix(input, &result);
  }

  result.finalize();

  for (const auto &vertex : this->mesh.localVertices()) {
    if (vertex.index() != vertexIndex) {
      continue;
    }

    const auto range = vertex.globalDofLineRange();
    auto reference = value_type{1.};

    for (auto i = range.first; i < range.second; ++i) {
      for (auto j = range.first; j < range.second; ++j) {
        EXPECT_THAT(result, AlmostEqIfLocal(i, j, reference));
        reference++;
      }
    }
  }
}

template <typename T> struct Mesh_IgnoreVertexTest : Test {
  using mesh_type = Mesh<T>;

  using size_type = typename mesh_type::size_type;
  using value_type = typename mesh_type::value_type;

  static constexpr size_type totalNumberOfElements = 2;
  static constexpr size_type totalNumberOfVertices = 2;
  static constexpr size_type dof = 0;
  static constexpr size_type dimension = 1;

  using Connectivity = std::vector<std::array<size_type, 2>>;
  Connectivity connectivity = {{mesh_type::IGNORE_VERTEX_INDEX, 1}, {0, 1}};

  mesh_type mesh = mesh_type::fromConnectivity(dimension, connectivity,
                                               totalNumberOfVertices, dof, dof);
};

using IgnoreVertex_TestCases =
    Types<SequentialComputePolicy, ParallelComputePolicy>;
TYPED_TEST_CASE(Mesh_IgnoreVertexTest, IgnoreVertex_TestCases);

TYPED_TEST(Mesh_IgnoreVertexTest, element_has_only_not_ignored_vertex) {
  for (const auto &element : this->mesh.localElements()) {
    if (element.index() == 0) {
      ASSERT_THAT(element.numberOfVertices(), Eq(1));
      EXPECT_THAT(element.vertexIndices(), ElementsAre(Eq(1)));
    }
  }
}

} // namespace
} // namespace cpppetsc
} // namespace ae108
