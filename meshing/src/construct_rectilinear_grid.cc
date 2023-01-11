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

#include "ae108/meshing/construct_rectilinear_grid.h"
#include <Eigen/Dense>

namespace ae108 {
namespace meshing {

template <>
std::tuple<std::vector<std::array<double, 1>>,
           std::vector<std::array<std::size_t, 2>>>
construct_rectilinear_grid(
    const std::array<std::vector<double>, 1> &instances) noexcept {

  const auto generate_points = [&]() {
    std::vector<std::array<double, 1>> points;
    points.reserve(instances[0].size());

    for (const auto &x : instances[0])
      points.push_back({{x}});

    return points;
  };

  const auto step_to_index = [&](const std::size_t step_x) { return step_x; };

  const auto generate_edges = [&]() {
    std::vector<std::array<std::size_t, 2>> edges;
    edges.reserve(instances[0].size() - 1);
    std::array<std::size_t, 1> step = {{0}};

    // edges in x
    for (step[0] = 0; step[0] < instances[0].size() - 1; ++step[0])
      edges.push_back({{step_to_index(step[0]), step_to_index(step[0] + 1)}});

    return edges;
  };

  return {generate_points(), generate_edges()};
}

template <>
std::tuple<std::vector<std::array<double, 2>>,
           std::vector<std::array<std::size_t, 2>>>
construct_rectilinear_grid(
    const std::array<std::vector<double>, 2> &instances) noexcept {

  const auto generate_points = [&]() {
    std::vector<std::array<double, 2>> points;
    points.reserve(instances[0].size() * instances[1].size());

    for (const auto &x : instances[0])
      for (const auto &y : instances[1])
        points.push_back({{x, y}});

    return points;
  };

  const auto step_to_index = [&](const std::size_t step_x,
                                 const std::size_t step_y) {
    return step_x * instances[1].size() + step_y;
  };

  const auto generate_edges = [&]() {
    std::vector<std::array<std::size_t, 2>> edges;
    edges.reserve((instances[0].size() - 1) * instances[1].size() +
                  instances[0].size() * (instances[1].size() - 1));

    std::array<std::size_t, 2> step = {{0, 0}};

    // edges in y
    for (step[0] = 0; step[0] < instances[0].size(); ++step[0])
      for (step[1] = 0; step[1] < instances[1].size() - 1; ++step[1])
        edges.push_back({{step_to_index(step[0], step[1]),
                          step_to_index(step[0], step[1] + 1)}});

    // edges in x
    for (step[0] = 0; step[0] < instances[0].size() - 1; ++step[0])
      for (step[1] = 0; step[1] < instances[1].size(); ++step[1])
        edges.push_back({{step_to_index(step[0], step[1]),
                          step_to_index(step[0] + 1, step[1])}});
    return edges;
  };

  return {generate_points(), generate_edges()};
}

template <>
std::tuple<std::vector<std::array<double, 3>>,
           std::vector<std::array<std::size_t, 2>>>
construct_rectilinear_grid(
    const std::array<std::vector<double>, 3> &instances) noexcept {

  const auto generate_points = [&]() {
    std::vector<std::array<double, 3>> points;
    points.reserve(instances[0].size() * instances[1].size() *
                   instances[2].size());

    for (const auto &x : instances[0])
      for (const auto &y : instances[1])
        for (const auto &z : instances[2])
          points.push_back({{x, y, z}});

    return points;
  };

  const auto step_to_index = [&](const std::size_t step_x,
                                 const std::size_t step_y,
                                 const std::size_t step_z) {
    return step_x * instances[1].size() * instances[2].size() +
           step_y * instances[2].size() + step_z;
  };

  const auto generate_edges = [&]() {
    std::vector<std::array<std::size_t, 2>> edges;
    edges.reserve((instances[0].size() - 1) * (instances[1].size() - 1) *
                      instances[2].size() +
                  instances[0].size() * (instances[1].size() - 1) *
                      (instances[2].size() - 1) +
                  (instances[0].size() - 1) * instances[1].size() *
                      (instances[2].size() - 1));

    std::array<std::size_t, 3> step = {{0, 0, 0}};

    // edges in z
    for (step[0] = 0; step[0] < instances[0].size(); ++step[0])
      for (step[1] = 0; step[1] < instances[1].size(); ++step[1])
        for (step[2] = 0; step[2] < instances[2].size() - 1; ++step[2])
          edges.push_back({{step_to_index(step[0], step[1], step[2]),
                            step_to_index(step[0], step[1], step[2] + 1)}});

    // edges in y
    for (step[0] = 0; step[0] < instances[0].size(); ++step[0])
      for (step[1] = 0; step[1] < instances[1].size() - 1; ++step[1])
        for (step[2] = 0; step[2] < instances[2].size(); ++step[2])
          edges.push_back({{step_to_index(step[0], step[1], step[2]),
                            step_to_index(step[0], step[1] + 1, step[2])}});

    // edges in x
    for (step[0] = 0; step[0] < instances[0].size() - 1; ++step[0])
      for (step[1] = 0; step[1] < instances[1].size(); ++step[1])
        for (step[2] = 0; step[2] < instances[2].size(); ++step[2])
          edges.push_back({{step_to_index(step[0], step[1], step[2]),
                            step_to_index(step[0] + 1, step[1], step[2])}});

    return edges;
  };

  return {generate_points(), generate_edges()};
}

} // namespace meshing
} // namespace ae108