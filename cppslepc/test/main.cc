// © 2021 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/cpppetsc/SequentialComputePolicy.h"
#include "ae108/cppslepc/Context.h"
#include <gmock/gmock.h>

struct Policy {
  static void handleError(const PetscErrorCode) {}
};

int main(int argc, char **argv) {
  const auto context = ae108::cppslepc::Context<Policy>{&argc, &argv};
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}