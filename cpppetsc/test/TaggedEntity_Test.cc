// © 2020 ETH Zurich, Mechanics and Materials Lab
// © 2020 California Institute of Technology
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

#include "ae108/cpppetsc/TaggedEntity.h"
#include <gmock/gmock.h>
#include <type_traits>

using testing::DoubleEq;
using testing::Return;
using testing::StrictMock;
using testing::Test;

namespace ae108 {
namespace cpppetsc {
namespace {

struct Callable_Mock {
  MOCK_CONST_METHOD1(call_const, double(double));
  MOCK_METHOD1(call, double(double));

  double operator()(double x) const { return call_const(x); }
  double operator()(double x) { return call(x); }
};

struct TaggedEntity_Test : Test {
  struct MyTag {};

  template <class T> double unwrapImplicitly(T t) { return t; }
};

TEST_F(TaggedEntity_Test, tag_type_is_the_tag_type) {
  static_assert(
      std::is_same<MyTag, TaggedEntity<double, MyTag>::tag_type>::value,
      "tag_type should be MyTag");
}

TEST_F(TaggedEntity_Test, value_type_is_the_entity_type) {
  static_assert(
      std::is_same<double, TaggedEntity<double, MyTag>::value_type>::value,
      "value_type should be double");
}

TEST_F(TaggedEntity_Test, value_type_is_the_entity_type_without_const) {
  static_assert(
      std::is_same<double,
                   TaggedEntity<const double, MyTag>::value_type>::value,
      "value_type should be double");
}

TEST_F(TaggedEntity_Test, value_type_is_the_entity_type_without_const_ref) {
  static_assert(
      std::is_same<double,
                   TaggedEntity<const double &, MyTag>::value_type>::value,
      "value_type should be double");
}

TEST_F(TaggedEntity_Test, constructing_arguments_are_forwarded) {
  const double value = 7.;
  TaggedEntity<double, MyTag> entity(value);
  EXPECT_THAT(entity.unwrap(), DoubleEq(value));
}

TEST_F(TaggedEntity_Test, unwrapping_const_value_works) {
  const double value = 7.;
  const TaggedEntity<double, MyTag> entity(value);
  EXPECT_THAT(entity.unwrap(), DoubleEq(value));
}

TEST_F(TaggedEntity_Test, unwrapping_rref_works) {
  const double value = 7.;
  EXPECT_THAT((TaggedEntity<double, MyTag>(value).unwrap()), DoubleEq(value));
}

TEST_F(TaggedEntity_Test, calling_function_works) {
  TaggedEntity<StrictMock<Callable_Mock>, MyTag> entity;

  const double value = 7.;
  const double returnValue = 3.;
  EXPECT_CALL(entity.unwrap(), call(value)).WillOnce(Return(returnValue));

  EXPECT_THAT(entity(value), DoubleEq(returnValue));
}

TEST_F(TaggedEntity_Test, calling_const_function_works) {
  const TaggedEntity<StrictMock<Callable_Mock>, MyTag> entity;

  const double value = 7.;
  const double returnValue = 3.;
  EXPECT_CALL(entity.unwrap(), call_const(value)).WillOnce(Return(returnValue));

  EXPECT_THAT(entity(value), DoubleEq(returnValue));
}

TEST_F(TaggedEntity_Test, implicit_conversion_to_const_ref_works) {
  const double value = 7.;
  const TaggedEntity<double, MyTag> entity(value);
  EXPECT_THAT(unwrapImplicitly<const double &>(entity), DoubleEq(value));
}

TEST_F(TaggedEntity_Test, implicit_conversion_to_ref_works) {
  const double value = 7.;
  TaggedEntity<double, MyTag> entity(value);
  EXPECT_THAT(unwrapImplicitly<double &>(entity), DoubleEq(value));
}

TEST_F(TaggedEntity_Test, implicit_conversion_to_rref_works) {
  const double value = 7.;
  EXPECT_THAT(unwrapImplicitly<double &&>(TaggedEntity<double, MyTag>(value)),
              DoubleEq(value));
}

TEST_F(TaggedEntity_Test, tagging_works) {
  const double value = 7.;
  const auto entity = tag<MyTag>(value);

  static_assert(
      std::is_same<double, TaggedEntity<double, MyTag>::value_type>::value,
      "value_type should be double");
  static_assert(
      std::is_same<MyTag, TaggedEntity<double, MyTag>::tag_type>::value,
      "tag_type should be MyTag");
  EXPECT_THAT(entity.unwrap(), DoubleEq(value));
}
} // namespace
} // namespace cpppetsc
} // namespace ae108
