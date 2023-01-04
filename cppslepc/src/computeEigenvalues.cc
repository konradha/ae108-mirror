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
computeEigenvalues(const cpppetsc::Matrix<cpppetsc::SequentialComputePolicy> &,
                   const cpppetsc::Matrix<cpppetsc::SequentialComputePolicy> &);

template std::vector<std::complex<typename LinearEigenvalueProblemSolver<
    cpppetsc::ParallelComputePolicy>::real_type>>
computeEigenvalues(const cpppetsc::Matrix<cpppetsc::ParallelComputePolicy> &A,
                   const cpppetsc::Matrix<cpppetsc::ParallelComputePolicy> &B);

} // namespace cppslepc
} // namespace ae108