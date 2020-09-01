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

#include "ae108/cpppetsc/IndexBasedAccessIterator.h"
#include <gmock/gmock.h>
#include <type_traits>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace cpppetsc {
namespace {

struct IndexBasedAccessIterator_Test : Test {
  std::vector<int> values = {1, 2, 3};

  using iterator = IndexBasedAccessIterator<std::vector<int>>;

  iterator begin{&values, 0};
  iterator end{&values, 3};
};

TEST_F(IndexBasedAccessIterator_Test, correct_typedefs) {
  static_assert(std::is_same<iterator::value_type, int>::value,
                "value_type is the range's value_type");
  static_assert(std::is_same<iterator::iterator_category,
                             std::bidirectional_iterator_tag>::value,
                "it's a bidirectional iterator");
  static_assert(std::is_same<iterator::difference_type, std::ptrdiff_t>::value,
                "difference_type is ptrdiff_t");
  static_assert(std::is_same<iterator::reference, int>::value,
                "references can bind to value_type");
}

TEST_F(IndexBasedAccessIterator_Test, dereferencing_works) {
  EXPECT_THAT(*begin, Eq(1));
}

TEST_F(IndexBasedAccessIterator_Test, assignment_works) {
  end = begin;
  EXPECT_THAT(end, Eq(begin));
}

TEST_F(IndexBasedAccessIterator_Test, prefix_increment_works) {
  const auto result = ++begin;
  EXPECT_THAT(*begin, Eq(2));
  EXPECT_THAT(*result, Eq(2));
}

TEST_F(IndexBasedAccessIterator_Test, postfix_increment_works) {
  const auto result = begin++;
  EXPECT_THAT(*begin, Eq(2));
  EXPECT_THAT(*result, Eq(1));
}

TEST_F(IndexBasedAccessIterator_Test, prefix_decrement_works) {
  ++begin;
  const auto result = --begin;
  EXPECT_THAT(*begin, Eq(1));
  EXPECT_THAT(*result, Eq(1));
}

TEST_F(IndexBasedAccessIterator_Test, postfix_decrement_works) {
  ++begin;
  const auto result = begin--;
  EXPECT_THAT(*begin, Eq(1));
  EXPECT_THAT(*result, Eq(2));
}

TEST_F(IndexBasedAccessIterator_Test, equal_comparison_works) {
  const auto copy = begin;
  EXPECT_THAT(begin == copy, Eq(true));
  EXPECT_THAT(begin == end, Eq(false));
}

TEST_F(IndexBasedAccessIterator_Test, not_equal_comparison_works) {
  const auto copy = begin;
  EXPECT_THAT(begin != copy, Eq(false));
  EXPECT_THAT(begin != end, Eq(true));
}

TEST_F(IndexBasedAccessIterator_Test, operator_arrow_works) {
  struct ABC {
    int m = 7;
  };
  std::vector<ABC> abc(1);
  IndexBasedAccessIterator<std::vector<ABC>> iter{&abc, 0};
  EXPECT_THAT(iter->m, Eq(7));
}
} // namespace
} // namespace cpppetsc
} // namespace ae108
