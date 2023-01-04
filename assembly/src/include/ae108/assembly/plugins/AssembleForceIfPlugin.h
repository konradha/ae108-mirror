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
#include <functional>

namespace ae108 {
namespace assembly {
namespace plugins {

/**
 * @brief Conditionally adds the force at a vertex to `output`.
 * The force is added if the `select` function returns true.
 */
DEFINE_CONST_ASSEMBLER_PLUGIN(
    AssembleForceIfPlugin, assembleForceIf,
    (const cpppetsc::local<vector_type> &displacements, const double time,
     std::function<bool(size_type elementIndex, size_type vertexIndex)> select,
     value_type *const output)) {
  typename element_type::NodalDisplacements elementInput;

  for (const auto &meshElement : this->assembler().meshElements()) {
    const auto &element = meshElement.meshView();
    utilities::groupElementDataPerVertex(element, displacements, &elementInput);

    const auto forces =
        meshElement.instance().computeForces(elementInput, time);

    auto forceIterator = forces.begin();
    for (auto &&vertex : element.vertexIndices()) {
      if (select(element.index(), vertex)) {
        for (auto index = decltype(forceIterator->size()){0};
             index < forceIterator->size(); ++index) {
          output[index] += (*forceIterator)[index];
        }
      }
      ++forceIterator;
    }
  }
}
} // namespace plugins
} // namespace assembly
} // namespace ae108
