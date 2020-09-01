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
#include <array>
#include <gmock/gmock.h>
#include <sstream>

using testing::DoubleEq;
using testing::Eq;
using testing::ExitedWithCode;
using testing::StrEq;
using testing::Test;

namespace ae108 {
namespace cmdline {
namespace {

struct CommandLineOptionParser_Test : Test {
  std::stringstream stream;
  CommandLineOptionParser parser{stream};
};

TEST_F(CommandLineOptionParser_Test, parses_double_argument) {
  constexpr int argc = 2;
  const std::array<const char *, argc> argv = {"test", "--value=.7"};

  double value = 0.;
  parser.withOption("value", &value).parse(argc, argv.data());

  EXPECT_THAT(value, DoubleEq(.7));
}

TEST_F(CommandLineOptionParser_Test, parses_short_options) {
  constexpr int argc = 3;
  const std::array<const char *, argc> argv = {"test", "-v", ".7"};

  double value = 0.;
  parser.withOption("value,v", &value).parse(argc, argv.data());

  EXPECT_THAT(value, DoubleEq(.7));
}

TEST_F(CommandLineOptionParser_Test, parses_bool_argument) {
  constexpr int argc = 2;
  const std::array<const char *, argc> argv = {"test", "--value=true"};

  bool value = false;
  parser.withOption("value", &value).parse(argc, argv.data());

  EXPECT_THAT(value, Eq(true));
}

TEST_F(CommandLineOptionParser_Test, parses_int_argument) {
  constexpr int argc = 2;
  const std::array<const char *, argc> argv = {"test", "--value=7"};

  int value = 0;
  parser.withOption("value", &value).parse(argc, argv.data());

  EXPECT_THAT(value, Eq(7));
}

TEST_F(CommandLineOptionParser_Test, does_not_change_value_if_missing) {
  constexpr int argc = 1;
  const std::array<const char *, argc> argv = {"test"};

  int value = 3;
  parser.withOption("value", &value).parse(argc, argv.data());

  EXPECT_THAT(value, Eq(3));
}

TEST_F(CommandLineOptionParser_Test, ignores_invalid_arguments) {
  constexpr int argc = 3;
  const std::array<const char *, argc> argv = {"test", "--unknown=123",
                                               "--known=1"};

  int known = 0;
  parser.withOption("known", &known).parse(argc, argv.data());
  EXPECT_THAT(known, Eq(1));
}

TEST_F(CommandLineOptionParser_Test, does_not_warn_for_no_unknown_args) {
  constexpr int argc = 2;
  const std::array<const char *, argc> argv = {"test", "--known=1"};

  int known = 0;
  parser.withOption("known", &known).parse(argc, argv.data());
  EXPECT_THAT(stream.str(), StrEq(""));
}

TEST_F(CommandLineOptionParser_Test, warns_about_single_unknown_argument) {
  constexpr int argc = 2;
  const std::array<const char *, argc> argv = {"test", "--unknown=123"};

  parser.parse(argc, argv.data());

  EXPECT_THAT(stream.str(), StrEq("Warning: The following options were not "
                                  "recognized by cmdline: '--unknown=123'.\n"));
}

TEST_F(CommandLineOptionParser_Test, warns_about_two_unknown_arguments) {
  constexpr int argc = 3;
  const std::array<const char *, argc> argv = {"test", "--unknown=123",
                                               "unknown_2=a"};

  parser.parse(argc, argv.data());

  EXPECT_THAT(
      stream.str(),
      StrEq("Warning: The following options were not "
            "recognized by cmdline: '--unknown=123', 'unknown_2=a'.\n"));
}

using CommandLineOptionParser_DeathTest = CommandLineOptionParser_Test;

TEST_F(CommandLineOptionParser_DeathTest, dies_if_help_requested) {
  constexpr int argc = 2;
  const std::array<const char *, argc> argv = {"test", "--help"};

  CommandLineOptionParser parser{std::cerr};
  EXPECT_EXIT(parser.parse(argc, argv.data()), ExitedWithCode(0),
              "Command line options");
}

TEST_F(CommandLineOptionParser_DeathTest, help_text_can_be_specified) {
  constexpr int argc = 2;
  const std::array<const char *, argc> argv = {"test", "--help"};

  CommandLineOptionParser parser{std::cerr};
  int value = 0;
  EXPECT_EXIT(
      parser.withOption("value", "help text", &value).parse(argc, argv.data()),
      ExitedWithCode(0), "help text");
}

TEST_F(CommandLineOptionParser_DeathTest, rejects_invalid_options) {
  constexpr int argc = 2;
  const std::array<const char *, argc> argv = {"test", "--value=abc"};

  CommandLineOptionParser parser{std::cerr};
  int value = 0;
  EXPECT_EXIT(parser.withOption("value", &value).parse(argc, argv.data()),
              ExitedWithCode(1), ".*'--value'.*");
}
} // namespace
} // namespace cmdline
} // namespace ae108
