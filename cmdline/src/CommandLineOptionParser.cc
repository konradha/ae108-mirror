// © 2020 ETH Zurich, Mechanics and Materials Lab
// © 2020 California Institute of Technology
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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
