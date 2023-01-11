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

namespace ae108 {
namespace meshing {
namespace cppgmsh {

template <std::size_t Dimension> struct Node {

  using Index = std::size_t;
  using Position = std::array<double, Dimension>;

  Index id;
  Position position;

  bool operator<(const Node<Dimension> &another) const {
    return id < another.id;
  }

  bool operator==(const Node<Dimension> &another) const {
    return id == another.id;
  }
};

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108