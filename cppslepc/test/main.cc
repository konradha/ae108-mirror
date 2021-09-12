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