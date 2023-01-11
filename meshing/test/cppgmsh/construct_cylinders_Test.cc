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

#include "ae108/meshing/cppgmsh/Context.h"
#include "ae108/meshing/cppgmsh/construct_cylinders.h"
#include <gmock/gmock.h>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct construct_cylinders_Test : Test {};

TEST_F(construct_cylinders_Test, constructs_empty) {
  const auto gmshContext = Context(0, 0);
  std::vector<std::array<double, 3>> nodes{{0, 0, 0}, {1, 0, 0}};
  std::vector<std::array<std::size_t, 2>> connectivity{};
  std::vector<double> radii{};
  bool capped = false;

  const auto result = construct_cylinders(nodes, connectivity, radii, capped);
  ASSERT_THAT(result.size(), Eq(0));
}

TEST_F(construct_cylinders_Test, constructs_single_cylinder) {
  const auto gmshContext = Context(0, 0);
  std::vector<std::array<double, 3>> nodes{{0, 0, 0}, {1, 0, 0}};
  std::vector<std::array<std::size_t, 2>> connectivity{{0, 1}};
  std::vector<double> radii{1.};
  bool capped = false;

  const auto result = construct_cylinders(nodes, connectivity, radii, capped);
  ASSERT_THAT(result.size(), Eq(1));
}

TEST_F(construct_cylinders_Test, constructs_two_cylinders) {
  const auto gmshContext = Context(0, 0);
  std::vector<std::array<double, 3>> nodes{
      {0, 0, 0}, {1, 0, 0}, {0, 0, 1}, {1, 0, 1}};
  std::vector<std::array<std::size_t, 2>> connectivity{{0, 1}, {2, 3}};
  std::vector<double> radii{0.1, 0.1};
  bool capped = false;

  const auto result = construct_cylinders(nodes, connectivity, radii, capped);
  ASSERT_THAT(result.size(), Eq(2));
}

TEST_F(construct_cylinders_Test, constructs_capped_cylinder) {
  const auto gmshContext = Context(0, 0);
  std::vector<std::array<double, 3>> nodes{{0, 0, 0}, {1, 0, 0}};
  std::vector<std::array<std::size_t, 2>> connectivity{{0, 1}};
  std::vector<double> radii{1.};
  bool capped = true;

  const auto result = construct_cylinders(nodes, connectivity, radii, capped);
  ASSERT_THAT(result.size(), Eq(1));
}
} // namespace cppgmsh
} // namespace meshing
} // namespace ae108