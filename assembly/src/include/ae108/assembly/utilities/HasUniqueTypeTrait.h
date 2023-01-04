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

#include "ae108/assembly/utilities/IsSingleUniqueType.h"

namespace ae108 {
namespace assembly {
namespace utilities {

/**
 * @brief This type trait compares "TypeTrait<Assembler>::type" for all types
 * "Assembler" in Assemblers. If they are all the same, then HasUniqueTypeTrait
 * contains a typedef "type" equal to "TypeTrait<Assembler>::type" and inherits
 * from std::true_type. If this is not the case, then HasUniqueTypeTrait does
 * not contain such a typedef and inherits from std::false_type.
 */
template <template <class> class TypeTrait, class... Assemblers>
struct HasUniqueTypeTrait
    : IsSingleUniqueType<typename TypeTrait<Assemblers>::type...> {};
} // namespace utilities
} // namespace assembly
} // namespace ae108
