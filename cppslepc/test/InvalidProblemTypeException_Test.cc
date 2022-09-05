// © 2022 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/cppslepc/InvalidProblemTypeException.h"
#include <gmock/gmock.h>

using testing::StrEq;
using testing::Test;

namespace ae108 {
namespace cppslepc {
namespace {

struct InvalidProblemTypeException_Test : Test {
  InvalidProblemTypeException exception;
};

TEST_F(InvalidProblemTypeException_Test, what_describes_exception) {
  EXPECT_THAT(
      exception.what(),
      StrEq(
          "The problem type does not match the number of operators provided."));
}

} // namespace
} // namespace cppslepc
} // namespace ae108