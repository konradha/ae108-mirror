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
#include "ae108/assembly/utilities/serialize.h"
#include "ae108/cpppetsc/TaggedVector.h"
#include <vector>

namespace ae108 {
namespace assembly {
namespace plugins {

DEFINE_CONST_ASSEMBLER_PLUGIN(
    AssembleStiffnessMatrixPlugin, assembleStiffnessMatrix,
    (const cpppetsc::local<vector_type> &displacements, const double time,
     matrix_type *const output)) {
  typename element_type::NodalDisplacements elementInput;
  std::vector<value_type> outputBuffer;

  for (const auto &meshElement : this->assembler().meshElements()) {
    utilities::groupElementDataPerVertex(meshElement.meshView(), displacements,
                                         &elementInput);

    const auto matrix =
        meshElement.instance().computeStiffnessMatrix(elementInput, time);
    outputBuffer.resize(matrix.rows() * matrix.cols());
    utilities::serialize(matrix, outputBuffer.begin());

    meshElement.meshView().addElementMatrix(outputBuffer, output);
  }
}
} // namespace plugins
} // namespace assembly
} // namespace ae108
