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

#include "ae108/elements/TimoshenkoBeamElementWithMass.h"
#include "ae108/elements/tensor/as_matrix_of_rows.h"
#include <Eigen/Geometry>

namespace ae108 {
namespace elements {
namespace {

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

/**
 * @brief Computes the consistent mass matrix of a reference beam with the given
 * properties. Since this matrix is diagonal, only the diagonal of the matrix is
 * stored.
 */
template <std::size_t Dimension_>
tensor::Tensor<double, Dimension_ *(Dimension_ + 1),
               Dimension_ *(Dimension_ + 1)>
consistent_mass_matrix(
    const TimoshenkoBeamProperties<double, Dimension_> &properties,
    const double length) noexcept;

// refer to Gan, Buntara S. (2018), "An Isogeometric Approach to Beam
// Structures", p.227, doi.org/10.1007/978-3-319-56493-7
template <>
tensor::Tensor<double, 6, 6>
consistent_mass_matrix<2>(const TimoshenkoBeamProperties<double, 2> &properties,
                          const double length) noexcept {
  const auto L = length;
  const auto A = properties.area;
  const auto rho = properties.density;
  const auto E = properties.young_modulus;
  const auto G = properties.shear_modulus;
  const auto I_z = properties.area_moment_z;
  const auto k_y = properties.shear_correction_factor_y;

  const auto phi_y = 12. * E * I_z * k_y / A / G / L / L;
  const auto _ = 0.;

  const auto CT = [&](const double phi) {
    return rho * A * L / (1 + phi) / (1 + phi);
  };

  const auto t1 = [&](const double phi) {
    return CT(phi) * (13. / 35 + 7. / 10 * phi + 1. / 3 * phi * phi);
  };

  const auto t2 = [&](const double phi) {
    return CT(phi) * (11. / 210 + 11. / 120 * phi + 1. / 24 * phi * phi) * L;
  };

  const auto t3 = [&](const double phi) {
    return CT(phi) * (1. / 105 + 1. / 60 * phi + 1. / 120 * phi * phi) * L * L;
  };

  const auto t4 = [&](const double phi) {
    return CT(phi) * (9. / 70 + 3. / 10 * phi + 1. / 6 * phi * phi);
  };

  const auto t5 = [&](const double phi) {
    return CT(phi) * (13. / 420 + 3. / 40 * phi + 1. / 24 * phi * phi) * L;
  };

  const auto t6 = [&](const double phi) {
    return CT(phi) * (1. / 140 + 1. / 60 * phi + 1. / 120 * phi * phi) * L * L;
  };

  const auto C_R = [&](const double I, const double phi) {
    return rho * I / (1 + phi) / (1 + phi) / L;
  };

  const auto r1 = [&](const double I, const double phi) {
    return C_R(I, phi) * 6. / 5;
  };
  const auto r2 = [&](const double I, const double phi) {
    return C_R(I, phi) * (1. / 10 - 1. / 2 * phi) * L;
  };
  const auto r3 = [&](const double I, const double phi) {
    return C_R(I, phi) * (2. / 15 + 1. / 6 * phi + 1. / 3 * phi * phi) * L * L;
  };
  const auto r4 = [&](const double I, const double phi) {
    return C_R(I, phi) * (1. / 30 + 1. / 6 * phi - 1. / 6 * phi * phi) * L * L;
  };

  const auto TY1 = t1(phi_y);
  const auto TY2 = t2(phi_y);
  const auto TY3 = t3(phi_y);
  const auto TY4 = t4(phi_y);
  const auto TY5 = t5(phi_y);
  const auto TY6 = t6(phi_y);
  const auto X = 1. / 3 * rho * A * L;

  // clang-format off
  tensor::Tensor<double, 6, 6> MT = {{
                                      {{   X,    _,    _,  X/2,    _,    _}},
                                      {{   _,  TY1,  TY2,    _,  TY4, -TY5}},
                                      {{   _,    _,  TY3,    _,  TY5, -TY6}},
                                      {{   _,    _,    _,    X,    _,    _}},
                                      {{   _,    _,    _,    _,  TY1, -TY2}},
                                      {{   _,    _,    _,    _,    _,  TY3}},
                                    }};
  // clang-format on

  const auto RY1 = r1(I_z, phi_y);
  const auto RY2 = r2(I_z, phi_y);
  const auto RY3 = r3(I_z, phi_y);
  const auto RY4 = r4(I_z, phi_y);

  // clang-format off
  tensor::Tensor<double, 6, 6> MR = {{
                                      {{   _,    _,    _,    _,    _,    _}},
                                      {{   _,  RY1,  RY2,    _, -RY1,  RY2}},
                                      {{   _,    _,  RY3,    _, -RY2, -RY4}},
                                      {{   _,    _,    _,    _,    _,    _}},
                                      {{   _,    _,    _,    _,  RY1, -RY2}},
                                      {{   _,    _,    _,    _,    _,  RY3}},
                                    }};
  // clang-format on

  return [&] {
    tensor::Tensor<double, 6, 6> out;
    tensor::as_matrix_of_rows(&out) =
        tensor::as_matrix_of_rows(&MT) + tensor::as_matrix_of_rows(&MR);
    return out;
  }();
}

// refer to Felippa et al (2015), "Mass Matrix Templates: General Description
// and 1D Examples", p.46, http://dx.doi.org/10.1007/s11831-014-9108-x
template <>
tensor::Tensor<double, 12, 12>
consistent_mass_matrix<3>(const TimoshenkoBeamProperties<double, 3> &properties,
                          const double length) noexcept {
  const auto L = length;
  const auto A = properties.area;
  const auto rho = properties.density;
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

  // translational inertia components
  const auto CT = [&](const double phi) {
    return rho * A * L / (1 + phi) / (1 + phi);
  };

  const auto t1 = [&](const double phi) {
    return CT(phi) * (13. / 35 + 7. / 10 * phi + 1. / 3 * phi * phi);
  };

  const auto t2 = [&](const double phi) {
    return CT(phi) * (11. / 210 + 11. / 120 * phi + 1. / 24 * phi * phi) * L;
  };

  const auto t3 = [&](const double phi) {
    return CT(phi) * (1. / 105 + 1. / 60 * phi + 1. / 120 * phi * phi) * L * L;
  };

  const auto t4 = [&](const double phi) {
    return CT(phi) * (9. / 70 + 3. / 10 * phi + 1. / 6 * phi * phi);
  };

  const auto t5 = [&](const double phi) {
    return CT(phi) * (13. / 420 + 3. / 40 * phi + 1. / 24 * phi * phi) * L;
  };

  const auto t6 = [&](const double phi) {
    return CT(phi) * (1. / 140 + 1. / 60 * phi + 1. / 120 * phi * phi) * L * L;
  };

  const auto C_R = [&](const double I, const double phi) {
    return rho * I / (1 + phi) / (1 + phi) / L;
  };

  const auto r1 = [&](const double I, const double phi) {
    return C_R(I, phi) * 6. / 5;
  };
  const auto r2 = [&](const double I, const double phi) {
    return C_R(I, phi) * (1. / 10 - 1. / 2 * phi) * L;
  };
  const auto r3 = [&](const double I, const double phi) {
    return C_R(I, phi) * (2. / 15 + 1. / 6 * phi + 1. / 3 * phi * phi) * L * L;
  };
  const auto r4 = [&](const double I, const double phi) {
    return C_R(I, phi) * (1. / 30 + 1. / 6 * phi - 1. / 6 * phi * phi) * L * L;
  };

  const auto TY1 = t1(phi_y);
  const auto TZ1 = t1(phi_z);
  const auto TY2 = t2(phi_y);
  const auto TZ2 = t2(phi_z);
  const auto TY3 = t3(phi_y);
  const auto TZ3 = t3(phi_z);
  const auto TY4 = t4(phi_y);
  const auto TZ4 = t4(phi_z);
  const auto TY5 = t5(phi_y);
  const auto TZ5 = t5(phi_z);
  const auto TY6 = t6(phi_y);
  const auto TZ6 = t6(phi_z);

  const auto X = 1. / 3. * rho * A * L;
  const auto S = J_x / 3. * rho * L;
  const auto _ = 0.;

  // clang-format off
  tensor::Tensor<double, 12, 12> M_CT = 
  {{
      {{  X,   _,   _,   _,   _,   _, X/2,   _,   _,   _,   _,   _}},
      {{  _, TY1,   _,   _,   _, TY2,   _, TY4,   _,   _,   _,-TY5}},
      {{  _,   _, TZ1,   _,-TZ2,   _,   _,   _, TZ4,   _, TZ5,   _}},
      {{  _,   _,   _,   S,   _,   _,   _,   _,   _, S/2,   _,   _}},
      {{  _,   _,   _,   _, TZ3,   _,   _,   _,-TZ5,   _,-TZ6,   _}},
      {{  _,   _,   _,   _,   _, TY3,   _, TY5,   _,   _,   _,-TY6}},
      {{  _,   _,   _,   _,   _,   _,   X,   _,   _,   _,   _,   _}},
      {{  _,   _,   _,   _,   _,   _,   _, TY1,   _,   _,   _,-TY2}},
      {{  _,   _,   _,   _,   _,   _,   _,   _, TZ1,   _, TZ2,   _}},
      {{  _,   _,   _,   _,   _,   _,   _,   _,   _,   S,   _,   _}},
      {{  _,   _,   _,   _,   _,   _,   _,   _,   _,   _, TZ3,   _}},
      {{  _,   _,   _,   _,   _,   _,   _,   _,   _,   _,   _, TY3}},
  }};
  // clang-format on

  // translational inertia components

  const auto RY1 = r1(I_z, phi_y);
  const auto RZ1 = r1(I_y, phi_z);
  const auto RY2 = r2(I_z, phi_y);
  const auto RZ2 = r2(I_y, phi_z);
  const auto RY3 = r3(I_z, phi_y);
  const auto RZ3 = r3(I_y, phi_z);
  const auto RY4 = r4(I_z, phi_y);
  const auto RZ4 = r4(I_y, phi_z);

  // clang-format off
  tensor::Tensor<double, 12, 12> M_CR =
  {{
      {{  _,   _,   _,   _,   _,   _,   _,   _,   _,   _,   _,   _}},
      {{  _, RY1,   _,   _,   _, RY2,   _,-RY1,   _,   _,   _, RY2}},
      {{  _,   _, RZ1,   _,-RZ2,   _,   _,   _,-RZ1,   _,-RZ2,   _}},
      {{  _,   _,   _,   _,   _,   _,   _,   _,   _,   _,   _,   _}},
      {{  _,   _,   _,   _, RZ3,   _,   _,   _, RZ2,   _,-RZ4,   _}},
      {{  _,   _,   _,   _,   _, RY3,   _,-RY2,   _,   _,   _,-RY4}},
      {{  _,   _,   _,   _,   _,   _,   _,   _,   _,   _,   _,   _}},
      {{  _,   _,   _,   _,   _,   _,   _, RY1,   _,   _,   _,-RY2}},
      {{  _,   _,   _,   _,   _,   _,   _,   _, RZ1,   _, RZ2,   _}},
      {{  _,   _,   _,   _,   _,   _,   _,   _,   _,   _,   _,   _}},
      {{  _,   _,   _,   _,   _,   _,   _,   _,   _,   _, RZ3,   _}},
      {{  _,   _,   _,   _,   _,   _,   _,   _,   _,   _,   _, RY3}},
  }};
  // clang-format on

  return [&] {
    tensor::Tensor<double, 12, 12> M_CMM;
    tensor::as_matrix_of_rows(&M_CMM) =
        tensor::as_matrix_of_rows(&M_CT) + tensor::as_matrix_of_rows(&M_CR);
    return M_CMM;
  }();
}

template <std::size_t Dimension_>
Eigen::Matrix<double, Dimension_ *(Dimension_ + 1),
              Dimension_ *(Dimension_ + 1), Eigen::RowMajor>
timoshenko_beam_consistent_mass_matrix(
    const tensor::Tensor<double, Dimension_> &axis,
    const TimoshenkoBeamProperties<double, Dimension_> &properties) noexcept {
  const auto reference = consistent_mass_matrix<Dimension_>(
      properties, tensor::as_vector(&axis).norm());
  // reference only contains values in the "upper" section of the matrix;
  // selfadjointView interprets reference as a symmetric matrix
  const auto reference_matrix = tensor::as_matrix_of_rows(&reference)
                                    .template selfadjointView<Eigen::Upper>();

  const auto rotation = rotation_matrix<Dimension_>(axis);
  const auto rotation_matrix = tensor::as_matrix_of_rows(&rotation);

  return rotation_matrix.transpose() * reference_matrix * rotation_matrix;
}

template Eigen::Matrix<double, 6, 6, Eigen::RowMajor>
timoshenko_beam_consistent_mass_matrix(
    const tensor::Tensor<double, 2> &axis,
    const TimoshenkoBeamProperties<double, 2> &properties) noexcept;

template Eigen::Matrix<double, 12, 12, Eigen::RowMajor>
timoshenko_beam_consistent_mass_matrix(
    const tensor::Tensor<double, 3> &axis,
    const TimoshenkoBeamProperties<double, 3> &properties) noexcept;

} // namespace elements
} // namespace ae108