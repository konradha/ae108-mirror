// © 2021 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/cpppetsc/ParallelComputePolicy_fwd.h"
#include "ae108/cpppetsc/SequentialComputePolicy_fwd.h"
#include <memory>
#include <petscsys.h>

namespace ae108 {
namespace cpppetsc {

template <class Policy> class Context {
public:
  /**
   * @brief Initializes PETSc. PETSc is finalized automatically when the
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
   * @brief Initializes PETSc.
   */
  static void initialize(int *const argc, char ***const argv,
                         const std::string &database, const std::string &help);

  /**
   * @brief Finalizes PETSc.
   */
  static void finalize(void *);

  using Token = std::unique_ptr<Context, decltype(&finalize)>;
  Token token_;
};

extern template class Context<SequentialComputePolicy>;
extern template class Context<ParallelComputePolicy>;

} // namespace cpppetsc
} // namespace ae108

namespace ae108 {
namespace cpppetsc {

template <class Policy>
Context<Policy>::Context(int *const argc, char ***const argv,
                         const std::string &database, const std::string &help)
    : token_((initialize(argc, argv, database, help), this), &finalize) {}

template <class Policy>
void Context<Policy>::initialize(int *const argc, char ***const argv,
                                 const std::string &database,
                                 const std::string &help) {
  Policy::handleError(
      PetscInitialize(argc, argv, database.empty() ? nullptr : database.c_str(),
                      help.empty() ? nullptr : help.c_str()));
}

template <class Policy> void Context<Policy>::finalize(void *) {
  Policy::handleError(PetscFinalize());
}

} // namespace cpppetsc
} // namespace ae108