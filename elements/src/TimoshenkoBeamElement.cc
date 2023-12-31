// © 2020, 2021 ETH Zurich, Mechanics and Materials Lab
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
stiffness_matrix(const TimoshenkoBeamProperties<double, Dimension_> &properties,
                 const double length) noexcept;

// refer to Cook et. al (2002), "Concepts and applications of Finite Element
// Analysis", 4th ed., p.27
template <>
tensor::Tensor<double, 12, 12>
stiffness_matrix<3>(const TimoshenkoBeamProperties<double, 3> &properties,
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
stiffness_matrix<2>(const TimoshenkoBeamProperties<double, 2> &properties,
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
    const TimoshenkoBeamProperties<double, Dimension_> &properties) noexcept {
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
    const TimoshenkoBeamProperties<double, 2> &properties) noexcept;

template Eigen::Matrix<double, 12, 12, Eigen::RowMajor>
timoshenko_beam_stiffness_matrix(
    const tensor::Tensor<double, 3> &axis,
    const TimoshenkoBeamProperties<double, 3> &properties) noexcept;

/**
 * @brief Computes the lumped mass matrix of a reference beam with the given
 * properties. Since this matrix is diagonal, only the diagonal of the matrix is
 * stored.
 */
template <std::size_t Dimension_>
tensor::Tensor<double, Dimension_ *(Dimension_ + 1), 1> lumped_mass_matrix(
    const TimoshenkoBeamProperties<double, Dimension_> &properties,
    const double density, const double length) noexcept;

// refer to Felippa et al (2015), "Mass Matrix Templates: General Description
// and 1D Examples", p.45, http://dx.doi.org/10.1007/s11831-014-9108-x
template <>
tensor::Tensor<double, 6, 1>
lumped_mass_matrix<2>(const TimoshenkoBeamProperties<double, 2> &properties,
                      const double density, const double length) noexcept {
  const auto mass = length * properties.area * density;
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
                      const double density, const double length) noexcept {
  const auto mass = length * properties.area * density;
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
    const TimoshenkoBeamProperties<double, Dimension_> &properties,
    const double density) noexcept {
  const auto reference = lumped_mass_matrix<Dimension_>(
      properties, density, tensor::as_vector(&axis).norm());
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
    const TimoshenkoBeamProperties<double, 2> &properties,
    const double density) noexcept;

template Eigen::Matrix<double, 12, 12, Eigen::RowMajor>
timoshenko_beam_lumped_mass_matrix(
    const tensor::Tensor<double, 3> &axis,
    const TimoshenkoBeamProperties<double, 3> &properties,
    const double density) noexcept;

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
    const double density, const double length) noexcept;

// refer to Felippa et al (2015), "Mass Matrix Templates: General Description
// and 1D Examples", Eq. 154, http://dx.doi.org/10.1007/s11831-014-9108-x
tensor::Tensor<double, 6> translational_inertia_terms(const double rho,
                                                      const double A,
                                                      const double L,
                                                      const double phi) {
  const auto C_T = rho * A * L / (1 + phi) / (1 + phi);
  return {C_T * (13. / 35 + 7. / 10 * phi + 1. / 3 * phi * phi),
          C_T * (11. / 210 + 11. / 120 * phi + 1. / 24 * phi * phi) * L,
          C_T * (1. / 105 + 1. / 60 * phi + 1. / 120 * phi * phi) * L * L,
          C_T * (9. / 70 + 3. / 10 * phi + 1. / 6 * phi * phi),
          C_T * (13. / 420 + 3. / 40 * phi + 1. / 24 * phi * phi) * L,
          C_T * (1. / 140 + 1. / 60 * phi + 1. / 120 * phi * phi) * L * L};
}

// refer to Felippa et al (2015), "Mass Matrix Templates: General Description
// and 1D Examples", Eq. 154, http://dx.doi.org/10.1007/s11831-014-9108-x
tensor::Tensor<double, 4> rotational_inertia_terms(const double rho,
                                                   const double I,
                                                   const double L,
                                                   const double phi) {
  const auto C_R = rho * I / (1 + phi) / (1 + phi) / L;
  return {C_R * 6. / 5, C_R * (1. / 10 - 1. / 2 * phi) * L,
          C_R * (2. / 15 + 1. / 6 * phi + 1. / 3 * phi * phi) * L * L,
          C_R * (1. / 30 + 1. / 6 * phi - 1. / 6 * phi * phi) * L * L};
}

// refer to Felippa et al (2015), "Mass Matrix Templates: General Description
// and 1D Examples", eq. 154, http://dx.doi.org/10.1007/s11831-014-9108-x
template <>
tensor::Tensor<double, 6, 6>
consistent_mass_matrix<2>(const TimoshenkoBeamProperties<double, 2> &properties,
                          const double density, const double length) noexcept {
  const auto L = length;
  const auto A = properties.area;
  const auto rho = density;
  const auto E = properties.young_modulus;
  const auto G = properties.shear_modulus;
  const auto I_z = properties.area_moment_z;
  const auto k_y = properties.shear_correction_factor_y;

  const auto phi_y = 12. * E * I_z * k_y / A / G / L / L;

  const auto X = 1. / 3 * rho * A * L;
  const auto TY = translational_inertia_terms(rho, A, L, phi_y);
  const auto RY = rotational_inertia_terms(rho, I_z, L, phi_y);
  const auto _ = 0.;

  // clang-format off
  tensor::Tensor<double, 6, 6> MT = {{
                                      {{     X,     _,     _,   X/2,     _,     _}},
                                      {{     _, TY[0], TY[1],     _, TY[3],-TY[4]}},
                                      {{     _,     _, TY[2],     _, TY[4],-TY[5]}},
                                      {{     _,     _,     _,     X,     _,     _}},
                                      {{     _,     _,     _,     _, TY[0],-TY[1]}},
                                      {{     _,     _,     _,     _,     _, TY[2]}},
                                    }};

  tensor::Tensor<double, 6, 6> MR = {{
                                      {{     _,     _,     _,     _,     _,     _}},
                                      {{     _, RY[0], RY[1],     _,-RY[0], RY[1]}},
                                      {{     _,     _, RY[2],     _,-RY[1],-RY[3]}},
                                      {{     _,     _,     _,     _,     _,     _}},
                                      {{     _,     _,     _,     _, RY[0],-RY[1]}},
                                      {{     _,     _,     _,     _,     _,  RY[2]}},
                                    }};
  // clang-format on

  return [&] {
    tensor::Tensor<double, 6, 6> out;
    tensor::as_matrix_of_rows(&out) =
        tensor::as_matrix_of_rows(&MT) + tensor::as_matrix_of_rows(&MR);
    return out;
  }();
}

template <>
tensor::Tensor<double, 12, 12>
consistent_mass_matrix<3>(const TimoshenkoBeamProperties<double, 3> &properties,
                          const double density, const double length) noexcept {
  const auto L = length;
  const auto A = properties.area;
  const auto rho = density;
  const auto E = properties.young_modulus;
  const auto G = properties.shear_modulus;
  const auto I_z = properties.area_moment_z;
  const auto k_y = properties.shear_correction_factor_y;
  const auto I_y = properties.area_moment_y;
  const auto J_x = properties.polar_moment_x;
  const auto k_z = properties.shear_correction_factor_z;

  const auto phi_y = 12. * E * I_z * k_y / A / G / L / L;
  const auto phi_z = 12. * E * I_y * k_z / A / G / L / L;

  const auto TY = translational_inertia_terms(rho, A, L, phi_y);
  const auto TZ = translational_inertia_terms(rho, A, L, phi_z);
  const auto RY = rotational_inertia_terms(rho, I_z, L, phi_y);
  const auto RZ = rotational_inertia_terms(rho, I_y, L, phi_z);
  const auto X = 1. / 3. * rho * A * L;
  const auto S = J_x / 3. * rho * L;
  const auto _ = 0.;

  // clang-format off
  tensor::Tensor<double, 12, 12> M_CT = 
  {{
      {{    X,     _,     _,     _,     _,     _,   X/2,     _,     _,     _,     _,     _}},
      {{    _, TY[0],     _,     _,     _, TY[1],     _, TY[3],     _,     _,     _,-TY[4]}},
      {{    _,     _, TZ[0],     _,-TZ[1],     _,     _,     _, TZ[3],     _, TZ[4],     _}},
      {{    _,     _,     _,     S,     _,     _,     _,     _,     _,   S/2,     _,     _}},
      {{    _,     _,     _,     _, TZ[2],     _,     _,     _,-TZ[4],     _,-TZ[5],     _}},
      {{    _,     _,     _,     _,     _, TY[2],     _, TY[4],     _,     _,     _,-TY[5]}},
      {{    _,     _,     _,     _,     _,     _,     X,     _,     _,     _,     _,     _}},
      {{    _,     _,     _,     _,     _,     _,     _, TY[0],     _,     _,     _,-TY[1]}},
      {{    _,     _,     _,     _,     _,     _,     _,     _, TZ[0],     _, TZ[1],     _}},
      {{    _,     _,     _,     _,     _,     _,     _,     _,     _,     S,     _,     _}},
      {{    _,     _,     _,     _,     _,     _,     _,     _,     _,     _, TZ[2],     _}},
      {{    _,     _,     _,     _,     _,     _,     _,     _,     _,     _,     _, TY[2]}},
  }};

  tensor::Tensor<double, 12, 12> M_CR =
  {{
      {{    _,     _,     _,     _,     _,     _,     _,     _,     _,     _,     _,     _}},
      {{    _, RY[0],     _,     _,     _, RY[1],     _,-RY[0],     _,     _,     _, RY[1]}},
      {{    _,     _, RZ[0],     _,-RZ[1],     _,     _,     _,-RZ[0],     _,-RZ[1],     _}},
      {{    _,     _,     _,     _,     _,     _,     _,     _,     _,     _,     _,     _}},
      {{    _,     _,     _,     _, RZ[2],     _,     _,     _, RZ[1],     _,-RZ[3],     _}},
      {{    _,     _,     _,     _,     _, RY[2],     _,-RY[1],     _,     _,     _,-RY[3]}},
      {{    _,     _,     _,     _,     _,     _,     _,     _,     _,     _,     _,     _}},
      {{    _,     _,     _,     _,     _,     _,     _, RY[0],     _,     _,     _,-RY[1]}},
      {{    _,     _,     _,     _,     _,     _,     _,     _, RZ[0],     _, RZ[1],     _}},
      {{    _,     _,     _,     _,     _,     _,     _,     _,     _,     _,     _,     _}},
      {{    _,     _,     _,     _,     _,     _,     _,     _,     _,     _, RZ[2],     _}},
      {{    _,     _,     _,     _,     _,     _,     _,     _,     _,     _,     _, RY[2]}},
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
    const TimoshenkoBeamProperties<double, Dimension_> &properties,
    const double density) noexcept {
  const auto reference = consistent_mass_matrix<Dimension_>(
      properties, density, tensor::as_vector(&axis).norm());
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
    const TimoshenkoBeamProperties<double, 2> &properties,
    const double density) noexcept;

template Eigen::Matrix<double, 12, 12, Eigen::RowMajor>
timoshenko_beam_consistent_mass_matrix(
    const tensor::Tensor<double, 3> &axis,
    const TimoshenkoBeamProperties<double, 3> &properties,
    const double density) noexcept;

} // namespace elements
} // namespace ae108