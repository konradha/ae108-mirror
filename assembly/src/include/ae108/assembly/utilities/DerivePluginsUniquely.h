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
#include "ae108/assembly/utilities/DeriveUniquely.h"

namespace ae108 {
namespace assembly {

/**
 * @tparam PluginList A FeaturePlugins-list of plugins.
 * @tparam Assembler The class that parametrizes the plugins.
 * @brief The struct inherits all the plugins once for the given assembler.
 */
template <class Assembler, class PluginList> struct DerivePluginsUniquely;

} // namespace assembly
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

namespace ae108 {
namespace assembly {

template <class Assembler, template <class> class... Plugins>
struct DerivePluginsUniquely<Assembler, assembly::FeaturePlugins<Plugins...>>
    : utilities::DeriveUniquely<Plugins<Assembler>...> {};

} // namespace assembly
} // namespace ae108
