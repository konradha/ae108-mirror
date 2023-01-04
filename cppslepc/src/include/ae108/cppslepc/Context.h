// © 2021 ETH Zurich, Mechanics and Materials Lab
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

#pragma once

#include "ae108/cpppetsc/ParallelComputePolicy_fwd.h"
#include "ae108/cpppetsc/SequentialComputePolicy_fwd.h"
#include <memory>
#include <slepcsys.h>

namespace ae108 {
namespace cppslepc {

template <class Policy> class Context {
public:
  /**
   * @brief Initializes SLEPc. SLEPc is finalized automatically when the
   * instance of Context is destroyed.
   *
   * @param argc A pointer to `argc` as provided to main().
   * @param argv A pointer to `argv` as provided to main().
   * @param database The path to a PETSC database file.
   * @param help Additional help to print.
   */
  explicit Context(int *const argc, char ***const argv,
                   const std::string &database = {},
                   const std::string &help = {});

private:
  /**
   * @brief Initializes SLEPc.
   */
  static void initialize(int *const argc, char ***const argv,
                         const std::string &database, const std::string &help);

  /**
   * @brief Finalizes SLEPc.
   */
  static void finalize(void *);

  using Token = std::unique_ptr<Context, decltype(&finalize)>;
  Token token_;
};

extern template class Context<cpppetsc::SequentialComputePolicy>;
extern template class Context<cpppetsc::ParallelComputePolicy>;

} // namespace cppslepc
} // namespace ae108

namespace ae108 {
namespace cppslepc {

template <class Policy>
Context<Policy>::Context(int *const argc, char ***const argv,
                         const std::string &database, const std::string &help)
    : token_((initialize(argc, argv, database, help), this), &finalize) {}

template <class Policy>
void Context<Policy>::initialize(int *const argc, char ***const argv,
                                 const std::string &database,
                                 const std::string &help) {
  Policy::handleError(
      SlepcInitialize(argc, argv, database.empty() ? nullptr : database.c_str(),
                      help.empty() ? nullptr : help.c_str()));
}

template <class Policy> void Context<Policy>::finalize(void *) {
  Policy::handleError(SlepcFinalize());
}

} // namespace cppslepc
} // namespace ae108