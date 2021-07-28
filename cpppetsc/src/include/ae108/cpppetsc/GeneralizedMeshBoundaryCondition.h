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

#include <vector>

namespace ae108 {
namespace cpppetsc {

/**
 * @brief Defines a boundary condition that sets the target value
 * (specified by a vertex view and a degree of freedom index) to
 * the sum of scaled source values (specified by a factor, a vertex index
 * and a degree of freedom index).
 */
template <class MeshType> struct GeneralizedMeshBoundaryCondition {
  using value_type = typename MeshType::value_type;
  using size_type = typename MeshType::size_type;

  struct Item {
    size_type vertex;
    size_type dof;
  };

  template <class Item> struct ScaledItem {
    value_type factor;
    Item item;
  };

  using target_type = Item;
  using sources_type = std::vector<ScaledItem<Item>>;

  target_type target;
  sources_type source;
  value_type offset;
};

} // namespace cpppetsc
} // namespace ae108