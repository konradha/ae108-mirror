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

#include "ae108/meshing/cppgmsh/as_affine_transform.h"
#include <Eigen/Dense>

namespace ae108 {
namespace meshing {
namespace cppgmsh {

std::array<double, 16>
as_affine_transform(const std::array<double, 3> &translation) noexcept {

  const auto transform =
      Eigen::Transform<double, 3, Eigen::Affine, Eigen::RowMajor>(
          Eigen::Translation3d(
              Eigen::Map<const Eigen::Vector3d>(translation.data())));

  return [transform]() {
    assert(transform.matrix().size() == 16);
    std::array<double, 16> out;
    std::move(transform.data(), transform.data() + 16, out.begin());
    return out;
  }();
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108