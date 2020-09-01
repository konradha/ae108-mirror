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
