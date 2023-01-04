// © 2020 ETH Zurich, Mechanics and Materials Lab
// © 2020 California Institute of Technology
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

#include "ae108/cppptest/isLocal.h"
#include <algorithm>
#include <cmath>
#include <complex>
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
inline bool almost_equal_to_reference(
    const double value, const double reference,
    ::testing::MatchResultListener *const stream) noexcept {
  return ::testing::ExplainMatchResult(::testing::DoubleEq(reference), value,
                                       stream);
}

/**
 * @brief Returns true if and only if value is approximately equal to the
 * reference.
 */
inline bool almost_equal_to_reference(
    const std::complex<double> &value, const std::complex<double> &reference,
    ::testing::MatchResultListener *const stream) noexcept {
  return almost_equal_to_reference(value.real(), reference.real(), stream) &&
         almost_equal_to_reference(value.imag(), reference.imag(), stream);
}

/**
 * @brief Returns true if and only if value is approximately equal to the
 * reference.
 */
inline bool
near_to_reference(const double value, const double reference,
                  const double tolerance,
                  ::testing::MatchResultListener *const stream) noexcept {
  return ::testing::ExplainMatchResult(
      ::testing::DoubleNear(reference, tolerance), value, stream);
}

/**
 * @brief Returns true if and only if value is approximately equal to the
 * reference.
 */
inline bool
near_to_reference(const std::complex<double> &value,
                  const std::complex<double> &reference, const double tolerance,
                  ::testing::MatchResultListener *const stream) noexcept {
  return near_to_reference(value.real(), reference.real(), tolerance, stream) &&
         near_to_reference(value.imag(), reference.imag(), tolerance, stream);
}

/**
 * @brief Check that the value at index row is almost equal to reference (if it
 * is available locally).
 */
MATCHER_P2(ScalarEqIfLocal, row, reference,
           std::string(negation ? "not " : "") + "equal to " +
               ::testing::PrintToString(reference) + " at (" +
               ::testing::PrintToString(row) + ")") {
  return !isLocal(arg, row) ||
         almost_equal_to_reference(arg(row), reference, result_listener);
}

/**
 * @brief Check that the value at index (row, col) is almost equal to reference
 * (if it is available locally).
 */
MATCHER_P3(ScalarEqIfLocal, row, col, reference,
           std::string(negation ? "not " : "") + "equal to " +
               ::testing::PrintToString(reference) + " at (" +
               ::testing::PrintToString(row) + ", " +
               ::testing::PrintToString(col) + ")") {
  return !isLocal(arg, row) ||
         almost_equal_to_reference(arg(row, col), reference, result_listener);
}

/**
 * @brief Check that the value at index row is almost equal to reference (if it
 * is available locally).
 */
MATCHER_P3(ScalarNearIfLocal, row, reference, tolerance,
           std::string(negation ? "not " : "") + "equal to " +
               ::testing::PrintToString(reference) + " at (" +
               ::testing::PrintToString(row) + ")") {
  return !isLocal(arg, row) ||
         near_to_reference(arg(row), reference, tolerance, result_listener);
}

/**
 * @brief Check that the value at index (row, col) is almost equal to reference
 * (if it is available locally).
 */
MATCHER_P4(ScalarNearIfLocal, row, col, reference, tolerance,
           std::string(negation ? "not " : "") + "equal to " +
               ::testing::PrintToString(reference) + " at (" +
               ::testing::PrintToString(row) + ", " +
               ::testing::PrintToString(col) + ")") {
  return !isLocal(arg, row) || near_to_reference(arg(row, col), reference,
                                                 tolerance, result_listener);
}

/**
 * @brief Check that the value has the given address.
 */
MATCHER_P(AddressEq, reference,
          "address " + std::string(negation ? "not " : "") + "equal to " +
              ::testing::PrintToString(reference)) {
  return ::testing::ExplainMatchResult(::testing::Eq(reference), &arg,
                                       result_listener);
}

/**
 * @brief Check that the value has almost the same values as reference.
 */
MATCHER_P(ScalarEq, reference,
          "value " + std::string(negation ? "not " : "") + "equal to " +
              ::testing::PrintToString(reference)) {
  return almost_equal_to_reference(arg, reference, result_listener);
}

/**
 * @brief Check that the value has almost the same values as reference.
 */
MATCHER(ScalarEq,
        "value " + std::string(negation ? "not " : "") + "equal to reference") {
  return almost_equal_to_reference(std::get<0>(arg), std::get<1>(arg),
                                   result_listener);
}

/**
 * @brief Check that the value has almost the same values as reference.
 */
MATCHER_P2(ScalarNear, reference, tolerance,
           "value " + std::string(negation ? "not " : "") + "equal to " +
               ::testing::PrintToString(reference)) {
  return near_to_reference(arg, reference, tolerance, result_listener);
}
} // namespace
} // namespace cppptest
} // namespace ae108
