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

#include "ae108/cppptest/isLocal.h"
#include <algorithm>
#include <cmath>
#include <gmock/gmock.h>
#include <limits>
#include <string>

namespace ae108 {
namespace cppptest {
namespace {

/**
 * @brief Returns true if and only if value is approximately equal to the
 * reference.
 */
inline bool almost_equal_to_reference(const double value,
                                      const double reference) noexcept {
  const auto tolerance = std::numeric_limits<double>::epsilon();
  return std::abs(value - reference) <=
         tolerance * std::max(tolerance, std::abs(reference));
}

/**
 * @brief Check that the value at index row is almost equal to reference (if it
 * is available locally).
 */
MATCHER_P2(AlmostEqIfLocal, row, reference,
           std::string(negation ? "not " : "") + "equal to " +
               ::testing::PrintToString(reference) + " at (" +
               ::testing::PrintToString(row) + ")") {
  return !isLocal(arg, row) || almost_equal_to_reference(arg(row), reference);
}

/**
 * @brief Check that the value at index (row, col) is almost equal to reference
 * (if it is available locally).
 */
MATCHER_P3(AlmostEqIfLocal, row, col, reference,
           std::string(negation ? "not " : "") + "equal to " +
               ::testing::PrintToString(reference) + " at (" +
               ::testing::PrintToString(row) + ", " +
               ::testing::PrintToString(col) + ")") {
  return !isLocal(arg, row) ||
         almost_equal_to_reference(arg(row, col), reference);
}

/**
 * @brief Check that the value has the given address.
 */
MATCHER_P(AddressEq, reference,
          "address " + std::string(negation ? "not " : "") + "equal to " +
              ::testing::PrintToString(reference)) {
  return reference == &arg;
}
} // namespace
} // namespace cppptest
} // namespace ae108