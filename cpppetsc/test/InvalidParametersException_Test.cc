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
