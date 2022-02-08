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

#include "ae108/cppslepc/computeEigenvalues.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"

namespace ae108 {
namespace cppslepc {

template std::vector<std::complex<typename LinearEigenvalueProblemSolver<
    cpppetsc::SequentialComputePolicy>::real_type>>
computeEigenvalues(const cpppetsc::Matrix<cpppetsc::SequentialComputePolicy> &);

template std::vector<std::complex<typename LinearEigenvalueProblemSolver<
    cpppetsc::ParallelComputePolicy>::real_type>>
computeEigenvalues(const cpppetsc::Matrix<cpppetsc::ParallelComputePolicy> &);

template std::vector<std::complex<typename LinearEigenvalueProblemSolver<
    cpppetsc::SequentialComputePolicy>::real_type>>
computeGeneralizedEigenvalues(
    const cpppetsc::Matrix<cpppetsc::SequentialComputePolicy> &,
    const cpppetsc::Matrix<cpppetsc::SequentialComputePolicy> &,
    const std::size_t number_of_eigenvalues);

template std::vector<std::complex<typename LinearEigenvalueProblemSolver<
    cpppetsc::ParallelComputePolicy>::real_type>>
computeGeneralizedEigenvalues(
    const cpppetsc::Matrix<cpppetsc::ParallelComputePolicy> &A,
    const cpppetsc::Matrix<cpppetsc::ParallelComputePolicy> &B,
    const std::size_t number_of_eigenvalues);

} // namespace cppslepc
} // namespace ae108