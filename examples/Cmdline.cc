// © 2020 ETH Zurich, Mechanics and Materials Lab
//
// This file is part of ae108.
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

#include <ae108/cmdline/CommandLineOptionParser.h>
#include <iostream>

int main(const int argc, const char *const *const argv) {
  namespace cmdline = ae108::cmdline;

  auto enable_greeting = bool{false};
  cmdline::CommandLineOptionParser(std::cerr)
      .withOption("enable_greeting,g", "Print a greeting.", &enable_greeting)
      .parse(argc, argv);

  if (enable_greeting) {
    std::cout << "Hello world!" << std::endl;
  }
}
