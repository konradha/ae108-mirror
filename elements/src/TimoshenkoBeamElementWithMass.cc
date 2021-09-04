// © 2020 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/elements/TimoshenkoBeamElementWithMass.h"
#include "ae108/elements/tensor/as_matrix_of_rows.h"
#include <Eigen/Geometry>

namespace ae108 {
namespace elements {
namespace {

/**
 * @brief Computes the lumped mass matrix of a reference beam with the given
 * properties. Since this matrix is diagonal, only the diagonal of the matrix is
 * stored.
 */
template <std::size_t Dimension_>
tensor::Tensor<double, Dimension_ *(Dimension_ + 1), 1> lumped_mass_matrix(
    const TimoshenkoBeamProperties<double, Dimension_> &properties,
    const double length) noexcept;

// refer to Felippa et al (2015), "Mass Matrix Templates: General Description
// and 1D Examples", p.45, http://dx.doi.org/10.1007/s11831-014-9108-x
template <>
tensor::Tensor<double, 6, 1>
lumped_mass_matrix<2>(const TimoshenkoBeamProperties<double, 2> &properties,
                      const double length) noexcept {
  const auto mass = length * properties.area * properties.density;
  const auto alpha =
      1. / 24; // ad hoc, refer to Cook et. al (2002), "Concepts and
               // applications of FiniteElement Analysis", 4th ed., p.378

  const auto A = mass / 2;
  const auto C = alpha * length * length * mass;

  return {{A, A, C, A, A, C}};
}

// refer to Felippa et al (2015), "Mass Matrix Templates: General Description
// and 1D Examples", p.45, http://dx.doi.org/10.1007/s11831-014-9108-x
template <>
tensor::Tensor<double, 12, 1>
lumped_mass_matrix<3>(const TimoshenkoBeamProperties<double, 3> &properties,
                      const double length) noexcept {
  const auto mass = length * properties.area * properties.density;
  const auto alpha =
      1. / 24; // ad hoc, refer to Cook et. al (2002), "Concepts and
               // applications of FiniteElement Analysis", 4th ed., p.378

  const auto A = mass / 2;
  const auto B = mass * properties.polar_moment_x / properties.area;
  const auto C = alpha * length * length * mass;

  return {{A, A, A, B, C, C, A, A, A, B, C, C}};
}

template <std::size_t Dimension_>
tensor::Tensor<double, Dimension_ *(Dimension_ + 1),
               Dimension_ *(Dimension_ + 1)>
rotation_matrix(const tensor::Tensor<double, Dimension_> &orientation) noexcept;

// refer to Cook et. al (2002), "Concepts and applications of Finite Element
// Analysis", 4th ed., p.32
template <>
tensor::Tensor<double, 12, 12>
rotation_matrix<3>(const tensor::Tensor<double, 3> &orientation) noexcept {
  // rotation that maps the normalized orientation vector to (1, 0, 0)
  const auto Lambda = Eigen::Quaternion<double>()
                          .FromTwoVectors(tensor::as_vector(&orientation),
                                          Eigen::Vector3d::UnitX())
                          .normalized()
                          .toRotationMatrix();

  auto result = tensor::Tensor<double, 12, 12>();
  auto result_matrix = tensor::as_matrix_of_rows(&result);

  result_matrix.block(0, 0, 3, 3) = Lambda;
  result_matrix.block(3, 3, 3, 3) = Lambda;
  result_matrix.block(6, 6, 3, 3) = Lambda;
  result_matrix.block(9, 9, 3, 3) = Lambda;

  return result;
}

// refer to Cook et. al (2002), "Concepts and applications of Finite Element
// Analysis", 4th ed., p.31
template <>
tensor::Tensor<double, 6, 6>
rotation_matrix<2>(const tensor::Tensor<double, 2> &orientation) noexcept {
  const auto normalized =
      tensor::as_vector(&orientation) / tensor::as_vector(&orientation).norm();

  // rotation that maps the normalized orientation vector to (1, 0)
  const tensor::Tensor<double, 2, 2> Lambda = {{
      {{normalized[0], normalized[1]}},
      {{-normalized[1], normalized[0]}},
  }};

  auto result = tensor::Tensor<double, 6, 6>();
  auto result_matrix = tensor::as_matrix_of_rows(&result);

  result_matrix.block(0, 0, 2, 2) = tensor::as_matrix_of_rows(&Lambda);
  result_matrix(2, 2) = 1.;
  result_matrix.block(3, 3, 2, 2) = tensor::as_matrix_of_rows(&Lambda);
  result_matrix(5, 5) = 1.;

  return result;
}
} // namespace

template <std::size_t Dimension_>
Eigen::Matrix<double, Dimension_ *(Dimension_ + 1),
              Dimension_ *(Dimension_ + 1), Eigen::RowMajor>
timoshenko_beam_lumped_mass_matrix(
    const tensor::Tensor<double, Dimension_> &axis,
    const TimoshenkoBeamProperties<double, Dimension_> &properties) noexcept {
  const auto reference = lumped_mass_matrix<Dimension_>(
      properties, tensor::as_vector(&axis).norm());
  // reference only contains values in the diagonal of the matrix;
  // asDiagonal interprets reference as a diagonal matrix
  const auto reference_matrix = tensor::as_vector(&reference).asDiagonal();

  const auto rotation = rotation_matrix<Dimension_>(axis);
  const auto rotation_matrix = tensor::as_matrix_of_rows(&rotation);

  return rotation_matrix.transpose() * reference_matrix * rotation_matrix;
}

template Eigen::Matrix<double, 6, 6, Eigen::RowMajor>
timoshenko_beam_lumped_mass_matrix(
    const tensor::Tensor<double, 2> &axis,
    const TimoshenkoBeamProperties<double, 2> &properties) noexcept;

template Eigen::Matrix<double, 12, 12, Eigen::RowMajor>
timoshenko_beam_lumped_mass_matrix(
    const tensor::Tensor<double, 3> &axis,
    const TimoshenkoBeamProperties<double, 3> &properties) noexcept;

} // namespace elements
} // namespace ae108