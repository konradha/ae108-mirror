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

#include "ae108/cpppetsc/LocalEntityIterator.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include <gmock/gmock.h>

namespace ae108 {
namespace cpppetsc {
namespace test {
namespace {

struct Mesh_Mock {
  using size_type = int;
  using value_type = double;
  using vector_type = std::vector<float>;
  using matrix_type = std::vector<std::vector<float>>;

  MOCK_CONST_METHOD1(vertices, std::vector<size_type>(size_type));
  MOCK_CONST_METHOD1(elementPointIndexToGlobalIndex, size_type(size_type));
  MOCK_CONST_METHOD1(vertexPointIndexToGlobalIndex, size_type(size_type));
  MOCK_CONST_METHOD1(numberOfVertices, size_type(size_type));

  MOCK_CONST_METHOD1(numberOfDofs, size_type(size_type));

  using Range = std::pair<size_type, size_type>;
  MOCK_CONST_METHOD1(localDofLineRange, Range(size_type));
  MOCK_CONST_METHOD1(globalDofLineRange, Range(size_type));

  MOCK_CONST_METHOD3(copyEntityData, void(size_type, const local<vector_type> &,
                                          std::vector<value_type> *));
  MOCK_CONST_METHOD3(addEntityData,
                     void(size_type, const std::vector<value_type> &,
                          local<vector_type> *));

  MOCK_CONST_METHOD3(addEntityMatrix,
                     void(size_type, const std::vector<value_type> &,
                          matrix_type *));
};

template <class Iterator> struct LocalEntityIterator_Test : ::testing::Test {
  Mesh_Mock mesh;

  typename Iterator::mesh_type::size_type begin_index = 1;
  Iterator begin{&mesh, begin_index};

  typename Iterator::mesh_type::size_type end_index = 7;
  Iterator end{&mesh, end_index};
};

TYPED_TEST_CASE_P(LocalEntityIterator_Test);

TYPED_TEST_P(LocalEntityIterator_Test, is_default_constructible) {
  EXPECT_NO_THROW(TypeParam());
}

TYPED_TEST_P(LocalEntityIterator_Test,
             eq_returns_false_for_different_iterators) {
  EXPECT_FALSE(this->begin == this->end);
}

TYPED_TEST_P(LocalEntityIterator_Test, eq_returns_true_for_same_iterator) {
  EXPECT_TRUE(this->begin == this->begin);
}

TYPED_TEST_P(LocalEntityIterator_Test,
             neq_returns_true_for_different_iterators) {
  EXPECT_TRUE(this->begin != this->end);
}

TYPED_TEST_P(LocalEntityIterator_Test, neq_returns_false_for_same_iterator) {
  EXPECT_FALSE(this->begin != this->begin);
}

TYPED_TEST_P(LocalEntityIterator_Test, is_copyconstrubtible) {
  TypeParam x(this->begin);
  EXPECT_THAT(x, ::testing::Eq(this->begin));
}

TYPED_TEST_P(LocalEntityIterator_Test, is_copyassignable) {
  TypeParam x;

  x = this->begin;

  EXPECT_THAT(x, ::testing::Eq(this->begin));
}

TYPED_TEST_P(LocalEntityIterator_Test, is_swappable) {
  using std::swap;

  TypeParam x(this->begin);
  TypeParam y(this->end);

  swap(x, y);

  EXPECT_THAT(x, ::testing::Eq(this->end));
  EXPECT_THAT(y, ::testing::Eq(this->begin));
}

TYPED_TEST_P(LocalEntityIterator_Test, is_dereferenceable) {
  EXPECT_NO_THROW(*this->begin);
}

TYPED_TEST_P(LocalEntityIterator_Test, is_preincrementable) {
  for (typename TypeParam::mesh_type::size_type i = this->begin_index;
       i < this->end_index; ++i)
    ++this->begin;

  EXPECT_THAT(this->begin, ::testing::Eq(this->end));
}

TYPED_TEST_P(LocalEntityIterator_Test, is_postincrementable) {
  for (typename TypeParam::mesh_type::size_type i = this->begin_index;
       i < this->end_index; ++i)
    this->begin++;

  EXPECT_THAT(this->begin, ::testing::Eq(this->end));
}

TYPED_TEST_P(LocalEntityIterator_Test, is_predecrementable) {
  for (typename TypeParam::mesh_type::size_type i = this->begin_index;
       i < this->end_index; ++i)
    --this->end;

  EXPECT_THAT(this->begin, ::testing::Eq(this->end));
}

TYPED_TEST_P(LocalEntityIterator_Test, is_postdecrementable) {
  for (typename TypeParam::mesh_type::size_type i = this->begin_index;
       i < this->end_index; ++i)
    this->end--;

  EXPECT_THAT(this->begin, ::testing::Eq(this->end));
}

TYPED_TEST_P(LocalEntityIterator_Test, preincrement_returns_new_iterator) {
  auto iter = this->begin;
  EXPECT_THAT(this->begin, ::testing::Not(::testing::Eq(++iter)));
}

TYPED_TEST_P(LocalEntityIterator_Test, postincrement_returns_old_iterator) {
  auto iter = this->begin;
  EXPECT_THAT(this->begin, ::testing::Eq(iter++));
}

TYPED_TEST_P(LocalEntityIterator_Test, predecrement_returns_new_iterator) {
  auto iter = this->end;
  EXPECT_THAT(this->end, ::testing::Not(::testing::Eq(--iter)));
}

TYPED_TEST_P(LocalEntityIterator_Test, postdecrement_returns_old_iterator) {
  auto iter = this->end;
  EXPECT_THAT(this->end, ::testing::Eq(iter--));
}

REGISTER_TYPED_TEST_CASE_P(
    LocalEntityIterator_Test, is_default_constructible,
    eq_returns_false_for_different_iterators, eq_returns_true_for_same_iterator,
    neq_returns_true_for_different_iterators,
    neq_returns_false_for_same_iterator, is_copyconstrubtible,
    is_copyassignable, is_swappable, is_dereferenceable, is_preincrementable,
    is_postincrementable, is_predecrementable, is_postdecrementable,
    preincrement_returns_new_iterator, postincrement_returns_old_iterator,
    predecrement_returns_new_iterator, postdecrement_returns_old_iterator);
} // namespace
} // namespace test
} // namespace cpppetsc
} // namespace ae108
