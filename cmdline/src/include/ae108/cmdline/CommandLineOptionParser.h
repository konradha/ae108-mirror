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

#pragma once

#include <boost/program_options.hpp>
#include <cassert>
#include <ostream>

namespace ae108 {
namespace cmdline {

class CommandLineOptionParser {
public:
  /**
   * @brief Create an instance of the class.
   * @param stream Messages (e.g. warnings) are printed to this stream.
   */
  explicit CommandLineOptionParser(std::ostream &stream);

  /**
   * @brief Add a parsed command line option to the parser. A "help" option does
   * not need to be provided since it is added automatically.
   *
   * @param name The name of the option, e.g. "verbose" (which will add an
   * option --verbose). Use "verbose,v" to also provide "-v".
   * @param help A description of the flag.
   * @param target The value of the flag (if any) will be written to this
   * variable. Must be a valid nonzero pointer.
   */
  template <typename T>
  CommandLineOptionParser &withOption(const char *name, const char *help,
                                      T *target);
  /**
   * @brief Calls withOption with no help text.
   */
  template <typename T>
  CommandLineOptionParser &withOption(const char *name, T *target);

  /**
   * @brief Parses the command line arguments. If "help" is requested, then help
   * is printed and the application exits with return code 0. If an error
   * occurs, then the error message is printed and the application exits with
   * return code 1.
   *
   * @remark Unknown options are ignored and a warning is printed.
   */
  void parse(int argc, const char *const *argv) const;

private:
  boost::program_options::options_description description_;
  std::ostream &stream_;
};
} // namespace cmdline
} // namespace ae108

/********************************************************************
 *  implementations
 *******************************************************************/

namespace ae108 {
namespace cmdline {

template <typename T>
CommandLineOptionParser &
CommandLineOptionParser::withOption(const char *const name,
                                    const char *const help, T *const target) {
  assert(target);

  description_.add_options()(name, ::boost::program_options::value<T>(target),
                             help);
  return *this;
}

template <typename T>
CommandLineOptionParser &
CommandLineOptionParser::withOption(const char *const name, T *const target) {
  return withOption(name, "", target);
}
} // namespace cmdline
} // namespace ae108
