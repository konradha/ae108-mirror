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

#include <memory>

namespace ae108 {
namespace meshing {
namespace cppgmsh {

class Context {
public:
  /**
   * @brief Initializes Gmsh. Gmsh is finalized automatically when the
   * instance of Context is destroyed.
   *
   * @param argc `argc` as provided to main().
   * @param argv `argv` as provided to main().
   * @param readConfigFiles Option to read config file.
   * @note https://gitlab.onelab.info/gmsh/gmsh/-/blob/gmsh_4_8_4/api/gmsh.h#L63
   */
  explicit Context(int const argc, char **const argv,
                   const bool readConfigFiles = true,
                   const bool verbose = false);

private:
  /**
   * @brief Initializes Gmsh.
   */
  static void initialize(int const argc, char **const argv,
                         const bool readConfigFiles, const bool verbose);

  /**
   * @brief Finalizes Gmsh.
   */
  static void finalize(void *);

  using Token = std::unique_ptr<Context, decltype(&finalize)>;
  Token token_;
};

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108