// © 2020 ETH Zurich, Mechanics and Materials Lab
// © 2020 California Institute of Technology
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

#include "ae108/cmdline/CommandLineOptionParser.h"
#include <cstdlib>
#include <string>
#include <vector>

namespace ae108 {
namespace cmdline {
CommandLineOptionParser::CommandLineOptionParser(std::ostream &stream)
    : description_("Command line options"), stream_(stream) {
  description_.add_options()("help,h", "Show this help.");
}

namespace {

void printUnrecognizedOptions(const std::vector<std::string> &options,
                              std::ostream &stream) {
  if (options.empty()) {
    return;
  }

  stream << "Warning: The following options were not recognized by cmdline: ";
  for (auto iterator = options.begin(); iterator != options.end(); ++iterator) {
    if (iterator != options.begin()) {
      stream << ", ";
    }
    stream << "'" << *iterator << "'";
  }
  stream << "." << std::endl;
}

} // namespace

void CommandLineOptionParser::parse(const int argc,
                                    const char *const *const argv) const {
  namespace options = ::boost::program_options;
  options::variables_map map;
  try {
    const auto parsedOptions = options::command_line_parser(argc, argv)
                                   .options(description_)
                                   .allow_unregistered()
                                   .run();

    printUnrecognizedOptions(
        options::collect_unrecognized(parsedOptions.options,
                                      options::include_positional),
        stream_);
    options::store(parsedOptions, map);
    options::notify(map);
  } catch (const options::error &exception) {
    stream_ << exception.what() << std::endl;
    std::exit(1);
  }
  if (map.count("help") > 0) {
    stream_ << description_ << std::flush;
    std::exit(0);
  }
}
} // namespace cmdline
} // namespace ae108
