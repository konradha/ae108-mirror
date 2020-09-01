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

#include "ae108/cpppetsc/LocalElementIterator.h"
#include "ae108/cpppetsc/test/LocalEntityIterator_Test.h"
#include "ae108/cppptest/Matchers.h"
#include <gmock/gmock.h>

using ae108::cppptest::AddressEq;
using testing::Eq;
using testing::Return;
using testing::Types;

namespace ae108 {
namespace cpppetsc {
namespace test {
namespace {

using IteratorsToTest = Types<LocalElementIterator<Mesh_Mock>>;
INSTANTIATE_TYPED_TEST_CASE_P(LocalElementIteratorTest,
                              LocalEntityIterator_Test, IteratorsToTest);
} // namespace
} // namespace test

namespace {

struct LocalElementIterator_Test
    : test::LocalEntityIterator_Test<LocalElementIterator<test::Mesh_Mock>> {
  using Iterator = LocalElementIterator<test::Mesh_Mock>;
};

TEST_F(LocalElementIterator_Test, accessing_index_works) {
  const auto value = Iterator::mesh_type::size_type{77};
  EXPECT_CALL(mesh, elementPointIndexToGlobalIndex(begin_index))
      .WillOnce(Return(value));

  EXPECT_THAT(begin->index(), Eq(value));
}

TEST_F(LocalElementIterator_Test, accessing_vertices_works) {
  const std::vector<Iterator::mesh_type::size_type> vertices = {3, 4, 5};
  EXPECT_CALL(mesh, vertices(begin_index)).WillOnce(Return(vertices));

  EXPECT_THAT(begin->vertexIndices(), Eq(vertices));
}

TEST_F(LocalElementIterator_Test, accessing_number_of_vertices_works) {
  const auto value = Iterator::mesh_type::size_type{77};
  EXPECT_CALL(mesh, numberOfVertices(begin_index)).WillOnce(Return(value));

  EXPECT_THAT(begin->numberOfVertices(), Eq(value));
}

TEST_F(LocalElementIterator_Test, dereferencing_works) {
  const auto value = Iterator::mesh_type::size_type{77};
  EXPECT_CALL(mesh, numberOfVertices(begin_index)).WillOnce(Return(value));

  EXPECT_THAT((*begin).numberOfVertices(), Eq(value));
}

TEST_F(LocalElementIterator_Test, copy_element_data_works) {
  TaggedEntity<test::Mesh_Mock::vector_type, LocalTag> from;
  std::vector<double> to;

  EXPECT_CALL(mesh, copyEntityData(begin_index, AddressEq(&from), &to))
      .Times(1);

  begin->copyElementData(from, &to);
}

TEST_F(LocalElementIterator_Test, add_element_data_works) {
  TaggedEntity<test::Mesh_Mock::vector_type, LocalTag> to;
  std::vector<double> from;

  EXPECT_CALL(mesh, addEntityData(begin_index, from, &to)).Times(1);

  begin->addElementData(from, &to);
}

TEST_F(LocalElementIterator_Test, add_element_matrix_works) {
  test::Mesh_Mock::matrix_type to;
  std::vector<double> from;

  EXPECT_CALL(mesh, addEntityMatrix(begin_index, from, &to)).Times(1);

  begin->addElementMatrix(from, &to);
}

TEST_F(LocalElementIterator_Test, number_of_dofs_works) {
  const test::Mesh_Mock::size_type dofs = 77;
  EXPECT_CALL(mesh, numberOfDofs(begin_index)).WillOnce(Return(dofs));

  EXPECT_THAT(begin->numberOfDofs(), Eq(dofs));
}

TEST_F(LocalElementIterator_Test, forClonedMesh_redirects_to_cloned_mesh) {
  const test::Mesh_Mock::size_type dofs = 77;

  test::Mesh_Mock newMesh;
  EXPECT_CALL(newMesh, numberOfDofs(begin_index)).WillOnce(Return(dofs));
  EXPECT_THAT(begin->forClonedMesh(&newMesh).numberOfDofs(), Eq(dofs));
}

TEST_F(LocalElementIterator_Test, local_dof_range_works) {
  const auto reference = std::make_pair(77, 777);

  EXPECT_CALL(mesh, localDofLineRange(begin_index)).WillOnce(Return(reference));
  EXPECT_THAT(begin->localDofLineRange(), Eq(reference));
}

TEST_F(LocalElementIterator_Test, global_dof_range_works) {
  const auto reference = std::make_pair(77, 777);

  EXPECT_CALL(mesh, globalDofLineRange(begin_index))
      .WillOnce(Return(reference));
  EXPECT_THAT(begin->globalDofLineRange(), Eq(reference));
}
} // namespace
} // namespace cpppetsc
} // namespace ae108
