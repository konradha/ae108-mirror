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

#include "ae108/cpppetsc/asSchurComplement.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/SequentialComputePolicy.h"

namespace ae108 {
namespace cpppetsc {

template Matrix<SequentialComputePolicy>
asSchurComplement(const Matrix<SequentialComputePolicy> *,
                  const Matrix<SequentialComputePolicy> *,
                  const Matrix<SequentialComputePolicy> *,
                  const Matrix<SequentialComputePolicy> *);

template Matrix<ParallelComputePolicy>
asSchurComplement(const Matrix<ParallelComputePolicy> *,
                  const Matrix<ParallelComputePolicy> *,
                  const Matrix<ParallelComputePolicy> *,
                  const Matrix<ParallelComputePolicy> *);

template Matrix<SequentialComputePolicy> asSchurComplement(
    const Matrix<SequentialComputePolicy> *,
    const std::vector<typename Matrix<SequentialComputePolicy>::size_type> &);

template Matrix<ParallelComputePolicy> asSchurComplement(
    const Matrix<ParallelComputePolicy> *,
    const std::vector<typename Matrix<ParallelComputePolicy>::size_type> &);

} // namespace cpppetsc
} // namespace ae108