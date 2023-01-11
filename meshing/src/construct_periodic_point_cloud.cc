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

#include "ae108/meshing/construct_periodic_point_cloud.h"
#include <Eigen/Dense>
#include <range/v3/view/cartesian_product.hpp>
#include <range/v3/view/zip.hpp>

namespace ae108 {
namespace meshing {

constexpr auto get_array = [](auto &&...x) {
  return std::array{std::forward<decltype(x)>(x)...};
};

template <std::size_t Dimension>
const auto construct_point_cloud =
    [](const std::array<std::array<double, Dimension>, Dimension> &translations,
       const std::array<double, Dimension> &origin, const auto &permutations) {
      std::vector<std::array<double, Dimension>> pointCloud(
          permutations.size());
      for (auto &&[point, tuple] :
           ranges::views::zip(pointCloud, permutations)) {
        const auto permutation = std::apply(get_array, tuple);
        Eigen::Map<Eigen::Matrix<double, Dimension, 1>>(point.data()) =
            Eigen::Map<const Eigen::Matrix<double, Dimension, 1>>(
                origin.data()) +
            (Eigen::Map<const Eigen::Matrix<double, Dimension, Dimension,
                                            Eigen::ColMajor>>(
                 translations.front().data()) *
             Eigen::Map<const Eigen::Matrix<double, Dimension, 1>>(
                 permutation.data()));
      }
      return pointCloud;
    };

template <>
std::vector<std::array<double, 3>> construct_periodic_point_cloud(
    const std::array<std::array<double, 3>, 3> &translations,
    std::array<double, 3> origin, std::vector<double> set) noexcept {
  return construct_point_cloud<3>(
      translations, origin, ranges::views::cartesian_product(set, set, set));
}

template <>
std::vector<std::array<double, 2>> construct_periodic_point_cloud(
    const std::array<std::array<double, 2>, 2> &translations,
    std::array<double, 2> origin, std::vector<double> set) noexcept {

  return construct_point_cloud<2>(translations, origin,
                                  ranges::views::cartesian_product(set, set));
}

template <>
std::vector<std::array<double, 1>> construct_periodic_point_cloud(
    const std::array<std::array<double, 1>, 1> &translations,
    std::array<double, 1> origin, std::vector<double> set) noexcept {
  return construct_point_cloud<1>(translations, origin,
                                  ranges::views::cartesian_product(set));
}

} // namespace meshing
} // namespace ae108