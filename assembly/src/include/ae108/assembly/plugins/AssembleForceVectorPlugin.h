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

#pragma once

#include "ae108/assembly/FeaturePlugin.h"
#include "ae108/assembly/utilities/groupElementDataPerVertex.h"
#include "ae108/cpppetsc/TaggedVector.h"

namespace ae108 {
namespace assembly {
namespace plugins {

DEFINE_CONST_ASSEMBLER_PLUGIN(
    AssembleForceVectorPlugin, assembleForceVector,
    (const cpppetsc::local<vector_type> &displacements, const double time,
     cpppetsc::local<vector_type> *const output)) {
  typename element_type::NodalDisplacements elementInput;

  for (const auto &meshElement : this->assembler().meshElements()) {
    utilities::groupElementDataPerVertex(meshElement.meshView(), displacements,
                                         &elementInput);

    const auto forces =
        meshElement.instance().computeForces(elementInput, time);

    utilities::ungroupElementDataPerVertex(meshElement.meshView(), forces,
                                           output);
  }
}
} // namespace plugins
} // namespace assembly
} // namespace ae108
