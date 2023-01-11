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

#include <array>
#include <vector>

namespace ae108 {
namespace meshing {

/**
 * @brief Returns a vector of points that is spanned by a set of discrete
 * translations around the origin. It is also known as Bravais lattice.
 *
 * In 3D, each point is defined by
 *
 *          R = n1 x a1 + n2 x a2 + n3 x a3,
 *
 * where ai is a translation vector and ni assumes any of the integers in set.
 *
 * @param translations An array containing the translation vectors.
 * @param origin Location of the origin.
 * @param set A set of integer permutations
 * @note https://en.wikipedia.org/wiki/Bravais_lattice
 */

template <std::size_t Dimension>
std::vector<std::array<double, Dimension>> construct_periodic_point_cloud(
    const std::array<std::array<double, Dimension>, Dimension> &translations,
    std::array<double, Dimension> origin,
    std::vector<double> set = {-1, 0, 1}) noexcept;

} // namespace meshing
} // namespace ae108