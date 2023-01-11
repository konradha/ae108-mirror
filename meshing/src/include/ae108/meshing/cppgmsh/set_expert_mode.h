// © 2022 ETH Zurich, Mechanics and Materials Lab
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

namespace ae108 {
namespace meshing {
namespace cppgmsh {

/**
 * @brief Enables/disables expert mode (to disable all the messages and user
 * input requests meant for inexperienced users)
 * @param enable Specifies if the expert mode should be enabled (true) or
 * disabled (false).
 */
void set_expert_mode(const bool enable = true) noexcept;

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108