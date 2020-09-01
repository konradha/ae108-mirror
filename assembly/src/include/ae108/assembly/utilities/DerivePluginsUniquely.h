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
