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

#include "ae108/elements/TimoshenkoBeamElement.h"
#include "ae108/elements/tensor/as_matrix_of_rows.h"
#include <Eigen/Geometry>

namespace ae108 {
namespace elements {
namespace {

/**
 * @brief Computes the stiffness matrix of a reference beam with the given
 * properties. Since this matrix is symmetric, only the "upper" section of the
 * matrix is filled with values.
 */
template <std::size_t Dimension_>
tensor::Tensor<double, Dimension_ *(Dimension_ + 1),
               Dimension_ *(Dimension_ + 1)>
stiffness_matrix(const Properties<double, Dimension_> &properties,
                 const double length) noexcept;

// refer to Cook et. al (2002), "Concepts and applications of Finite Element
// Analysis", 4th ed., p.27
template <>
tensor::Tensor<double, 12, 12>
stiffness_matrix<3>(const Properties<double, 3> &properties,
                    const double length) noexcept {
  const auto L = length;
  const auto A = properties.area;
  const auto E = properties.young_modulus;
  const auto G = properties.shear_modulus;
  const auto I_z = properties.area_moment_z;
  const auto k_y = properties.shear_correction_factor_y;
  const auto I_y = properties.area_moment_y;
  const auto J_x = properties.polar_moment_x;
  const auto k_z = properties.shear_correction_factor_z;

  const auto phi = [&](const double I, const double k) {
    return 12. * E * I * k / A / G / L / L;
  };
  const auto phi_y = phi(I_z, k_y);
  const auto phi_z = phi(I_y, k_z);

  const auto f1 = [&](const double I, const double phi) {
    return 12. * E * I / (1. + phi) / L / L / L;
  };
  const auto Y1 = f1(I_z, phi_y);
  const auto Z1 = f1(I_y, phi_z);

  const auto f2 = [&](const double I, const double phi) {
    return 6. * E * I / (1. + phi) / L / L;
  };
  const auto Y2 = f2(I_z, phi_y);
  const auto Z2 = f2(I_y, phi_z);

  const auto f3 = [&](const double I, const double phi) {
    return (4. + phi) * E * I / (1. + phi) / L;
  };
  const auto Y3 = f3(I_z, phi_y);
  const auto Z3 = f3(I_y, phi_z);

  const auto f4 = [&](const double I, const double phi) {
    return (2. - phi) * E * I / (1. + phi) / L;
  };
  const auto Y4 = f4(I_z, phi_y);
  const auto Z4 = f4(I_y, phi_z);

  const auto X = A * E / L;
  const auto S = G * J_x / L;
  const auto _ = 0.;

  // clang-format off
  return {{
      {{  X,   _,   _,   _,   _,   _,  -X,   _,   _,   _,   _,   _}},
      {{  _,  Y1,   _,   _,   _,  Y2,   _, -Y1,   _,   _,   _,  Y2}},
      {{  _,   _,  Z1,   _, -Z2,   _,   _,   _, -Z1,   _, -Z2,   _}},
      {{  _,   _,   _,   S,   _,   _,   _,   _,   _,  -S,   _,   _}},
      {{  _,   _,   _,   _,  Z3,   _,   _,   _,  Z2,   _,  Z4,   _}},
      {{  _,   _,   _,   _,   _,  Y3,   _, -Y2,   _,   _,   _,  Y4}},
      {{  _,   _,   _,   _,   _,   _,   X,   _,   _,   _,   _,   _}},
      {{  _,   _,   _,   _,   _,   _,   _,  Y1,   _,   _,   _, -Y2}},
      {{  _,   _,   _,   _,   _,   _,   _,   _,  Z1,   _,  Z2,   _}},
      {{  _,   _,   _,   _,   _,   _,   _,   _,   _,   S,   _,   _}},
      {{  _,   _,   _,   _,   _,   _,   _,   _,   _,   _,  Z3,   _}},
      {{  _,   _,   _,   _,   _,   _,   _,   _,   _,   _,   _,  Y3}},
  }};
  // clang-format on
}

// refer to Cook et. al (2002), "Concepts and applications of Finite Element
// Analysis", 4th ed., p.26
template <>
tensor::Tensor<double, 6, 6>
stiffness_matrix<2>(const Properties<double, 2> &properties,
                    const double length) noexcept {
  const auto L = length;
  const auto A = properties.area;
  const auto E = properties.young_modulus;
  const auto G = properties.shear_modulus;
  const auto I_z = properties.area_moment_z;
  const auto k_y = properties.shear_correction_factor_y;

  const auto phi_y = 12. * E * I_z * k_y / A / G / L / L;

  const auto Y1 = 12. * E * I_z / (1. + phi_y) / L / L / L;
  const auto Y2 = 6. * E * I_z / (1. + phi_y) / L / L;
  const auto Y3 = (4. + phi_y) * E * I_z / (1. + phi_y) / L;
  const auto Y4 = (2. - phi_y) * E * I_z / (1. + phi_y) / L;

  const auto X = A * E / L;
  const auto _ = 0.;

  // clang-format off
  return {{
      {{  X,   _,   _,  -X,   _,   _}},
      {{  _,  Y1,  Y2,   _, -Y1,  Y2}},
      {{  _,   _,  Y3,   _, -Y2,  Y4}},
      {{  _,   _,   _,   X,   _,   _}},
      {{  _,   _,   _,   _,  Y1, -Y2}},
      {{  _,   _,   _,   _,   _,  Y3}},
  }};
  // clang-format on
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
timoshenko_beam_stiffness_matrix(
    const tensor::Tensor<double, Dimension_> &axis,
    const Properties<double, Dimension_> &properties) noexcept {
  const auto reference =
      stiffness_matrix<Dimension_>(properties, tensor::as_vector(&axis).norm());
  // reference only contains values in the "upper" section of the matrix;
  // selfadjointView interprets reference as a symmetric matrix
  const auto reference_matrix = tensor::as_matrix_of_rows(&reference)
                                    .template selfadjointView<Eigen::Upper>();

  const auto rotation = rotation_matrix<Dimension_>(axis);
  const auto rotation_matrix = tensor::as_matrix_of_rows(&rotation);

  return rotation_matrix.transpose() * reference_matrix * rotation_matrix;
}

template Eigen::Matrix<double, 6, 6, Eigen::RowMajor>
timoshenko_beam_stiffness_matrix(
    const tensor::Tensor<double, 2> &axis,
    const Properties<double, 2> &properties) noexcept;

template Eigen::Matrix<double, 12, 12, Eigen::RowMajor>
timoshenko_beam_stiffness_matrix(
    const tensor::Tensor<double, 3> &axis,
    const Properties<double, 3> &properties) noexcept;

} // namespace elements
} // namespace ae108