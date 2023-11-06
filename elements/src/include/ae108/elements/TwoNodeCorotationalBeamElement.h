// © 2023 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/elements/ComputeEnergyTrait.h"
#include "ae108/elements/ComputeForcesTrait.h"
#include "ae108/elements/ComputeStiffnessMatrixTrait.h"
#include "ae108/elements/ElementBase.h"
#include "ae108/elements/tensor/as_vector.h"

#include "ae108/elements/TimoshenkoBeamElement.h"
// #include "ae108/elements/tensor/as_matrix_of_rows.h"


namespace ae108 {
namespace elements {

template <class RealType_, std::size_t Dimension_>
struct TwoNodeCorotationalBeamProperties;

template <class RealType_>
struct TwoNodeCorotationalBeamProperties<RealType_, 2> {
  using real_type = RealType_;
};



}
}