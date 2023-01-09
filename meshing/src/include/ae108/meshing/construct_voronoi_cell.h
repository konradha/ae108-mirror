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

/**
 * @brief Returns the boundary representation of the Voronoi cell spanned by a
 * point cloud around the origin. The Voronoi cell is defined as the locus of
 * all points closer to the origin than to any other point in the point cloud.
 *
 * @param point_cloud A vector of points, a.k.a. sites, forming the point cloud.
 * @param periodic_faces A vector containing all periodic face pairs.
 * @note https://en.wikipedia.org/wiki/Voronoi_diagram
 */

template <class SizeType, class ValueType, SizeType Dimension>
BoundaryRepresentation<SizeType, ValueType, Dimension> construct_voronoi_cell(
    const std::vector<std::array<ValueType, Dimension>> &point_cloud,
    std::vector<std::pair<SizeType, SizeType>> *periodic_faces =
        nullptr) noexcept;

} // namespace meshing
} // namespace ae108