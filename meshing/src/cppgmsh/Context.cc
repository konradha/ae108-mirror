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
#include <gmsh.h>

namespace ae108 {
namespace meshing {
namespace cppgmsh {

Context::Context(int const argc, char **const argv, const bool readConfigFiles,
                 const bool verbose)
    : token_((initialize(argc, argv, readConfigFiles, verbose), this),
             &finalize) {}

void Context::initialize(int const argc, char **const argv,
                         const bool readConfigFiles, const bool verbose) {
  gmsh::initialize(argc, argv, readConfigFiles);
  if (verbose)
    gmsh::option::setNumber("General.Verbosity", 5);
  else
    gmsh::option::setNumber("General.Verbosity", 0);
}

void Context::finalize(void *) { gmsh::finalize(); }

} // namespace cppgmsh
} // namespace meshing
} // namespace ae108