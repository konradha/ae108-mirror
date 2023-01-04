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

#pragma once

#include "ae108/assembly/FeaturePlugin.h"
#include "ae108/assembly/utilities/groupElementDataPerVertex.h"
#include "ae108/cpppetsc/TaggedVector.h"

namespace ae108 {
namespace assembly {
namespace plugins {

DEFINE_ASSEMBLER_PLUGIN(UpdateInternalVariablesPlugin, updateInternalVariables,
                        (const cpppetsc::local<vector_type> &displacements,
                         const double time)) {
  typename element_type::NodalDisplacements elementInput;
  for (auto &meshElement : this->assembler().meshElements()) {
    utilities::groupElementDataPerVertex(meshElement.meshView(), displacements,
                                         &elementInput);

    meshElement.instance().updateInternalVariables(elementInput, time);
  }
}
} // namespace plugins
} // namespace assembly
} // namespace ae108
