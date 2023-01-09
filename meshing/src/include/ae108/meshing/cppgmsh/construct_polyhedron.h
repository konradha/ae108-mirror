// © 2021 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/meshing/BoundaryRepresentation.h"

namespace ae108 {
namespace meshing {
namespace cppgmsh {

/**
 * @brief Constructs a solid polyhedron from a BoundaryRepresentation using
 * gmsh. Returns the gmsh entity.
 *
 * @param brep A boundary representation containing the polyhedron.
 *
 */
template <class SizeType, class ValueType>
std::pair<int, int> construct_polyhedron(
    const BoundaryRepresentation<SizeType, ValueType, 3> &brep) noexcept;

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108