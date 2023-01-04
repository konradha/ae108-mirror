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

#include "ae108/assembly/FeaturePlugins.h"
#include "ae108/assembly/plugins/AssembleEnergyPlugin.h"
#include "ae108/assembly/plugins/AssembleForceVectorPlugin.h"
#include "ae108/assembly/plugins/AssembleStiffnessMatrixPlugin.h"
#include "ae108/assembly/plugins/UpdateInternalVariablesPlugin.h"

namespace ae108 {
namespace assembly {

/**
 * @brief A list of default assembler plugins.
 */
using DefaultFeaturePlugins =
    FeaturePlugins<plugins::AssembleEnergyPlugin,
                   plugins::AssembleForceVectorPlugin,
                   plugins::AssembleStiffnessMatrixPlugin,
                   plugins::UpdateInternalVariablesPlugin>;
} // namespace assembly
} // namespace ae108
