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