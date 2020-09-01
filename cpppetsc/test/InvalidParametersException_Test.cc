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

#include "ae108/cpppetsc/InvalidParametersException.h"
#include <gmock/gmock.h>

using testing::StrEq;
using testing::Test;

namespace ae108 {
namespace cpppetsc {
namespace {

struct InvalidParametersException_Test : Test {
  InvalidParametersException exception;
};

TEST_F(InvalidParametersException_Test,
       what_states_that_parameters_are_invalid) {
  EXPECT_THAT(exception.what(),
              StrEq("The parameters that were provided are invalid."));
}

} // namespace
} // namespace cpppetsc
} // namespace ae108
