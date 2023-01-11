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

#include "ae108/meshing/cppgmsh/Node.h"
#include <gmock/gmock.h>

using testing::Eq;
using testing::Test;

namespace ae108 {
namespace meshing {
namespace cppgmsh {

struct Node_Test : Test {};

TEST_F(Node_Test, returns_self_equality) {
  const auto node = Node<3>{1, {0, 1, 2}};
  ASSERT_THAT(node == node, Eq(true));
}

TEST_F(Node_Test, returns_self_equality_2) {
  const auto node = Node<3>{1, {0, 1, 2}};
  ASSERT_THAT(node < node, Eq(false));
}

TEST_F(Node_Test, returns_correct_inequality_1) {
  const auto node1 = Node<3>{1, {0, 1, 2}};
  const auto node2 = Node<3>{2, {0, 1, 2}};
  ASSERT_THAT(node1 == node2, Eq(false));
}

TEST_F(Node_Test, returns_correct_inequality_2) {
  const auto node1 = Node<3>{1, {0, 1, 2}};
  const auto node2 = Node<3>{2, {0, 1, 2}};
  ASSERT_THAT(node1 < node2, Eq(true));
}

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108