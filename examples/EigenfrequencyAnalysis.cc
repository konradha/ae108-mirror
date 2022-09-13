// © 2020, 2021 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/assembly/Assembler.h"
#include "ae108/assembly/plugins/AssembleMassMatrixPlugin.h"
#include "ae108/cpppetsc/Mesh.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cpppetsc/computeElementsOfMatrix.h"
#include "ae108/cppslepc/Context.h"
#include "ae108/cppslepc/computeSmallestEigenvalues.h"
#include "ae108/elements/ElementWithMass.h"
#include "ae108/elements/TimoshenkoBeamElement.h"
#include "ae108/elements/compute_mass_matrix.h"
#include "ae108/elements/mesh/refine_segment_mesh.h"
#include "ae108/elements/tensor/as_vector.h"
#include <array>
#include <range/v3/view/drop.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/take.hpp>

using namespace ae108;
namespace rv = ranges::cpp20::views;

using Policy = cpppetsc::ParallelComputePolicy;
using Context = cppslepc::Context<Policy>;
using Mesh = cpppetsc::Mesh<Policy>;
using Vector = cpppetsc::Vector<Policy>;

// In this example we will calculate the eigenfrequencies of an unsupported
// Euler-Bernoulli beam and compare to the analytical solution.
// https://en.wikipedia.org/wiki/Euler-Bernoulli_beam_theory#Example:_unsupported_(free-free)_beam

// A beam with length 20 is discretized into 10 elements
//
// 0--2--3--4--5--6--7--8--9--10--1

constexpr auto coordinate_dimension = Mesh::size_type{2};

using Point = std::array<Mesh::real_type, coordinate_dimension>;

using TimoshenkoBeamElement =
    elements::TimoshenkoBeamElement<coordinate_dimension, Vector::value_type,
                                    Vector::real_type>;
using Element = elements::ElementWithMass<TimoshenkoBeamElement>;

using Properties =
    elements::TimoshenkoBeamProperties<Mesh::real_type, coordinate_dimension>;

// Let's not forget about the AssembleMassMatrixPlugin as we will need it to
// form the eigenvalue problem
using Plugins =
    assembly::FeaturePlugins<assembly::plugins::AssembleEnergyPlugin,
                             assembly::plugins::AssembleForceVectorPlugin,
                             assembly::plugins::AssembleStiffnessMatrixPlugin,
                             assembly::plugins::AssembleMassMatrixPlugin>;

using Assembler = assembly::Assembler<Element, Plugins, Policy>;

// We assume Euler-Bernoulli beam theory and hence chose the shear correction
// factor to be zero, such that the Timoshenko beam description reduces to the
// Euler-Bernoulli beam theory
constexpr Mesh::real_type young_modulus = 1.;
constexpr Mesh::real_type poisson_ratio = 0.3;
constexpr Mesh::real_type shear_modulus =
    young_modulus / (2 * (1 + poisson_ratio));
constexpr Mesh::real_type shear_correction_factor_y = 0.;
constexpr Mesh::real_type beam_length = 20.;
constexpr Mesh::real_type width = 1.;
constexpr Mesh::real_type thickness = 1.;
constexpr Mesh::real_type area = width * thickness;
constexpr Mesh::real_type density = 1.;
constexpr Mesh::real_type area_moment_z =
    thickness * width * width * width / 12.;

// It can be shown that the natural frequencies of an unsupported (free-free)
// Euler-Bernoulli beam are as follows:
const auto analytic_result = [](std::size_t n) {
  return std::pow((2 * n + 1) * M_PI / 2 / beam_length, 2) *
         std::sqrt(young_modulus * area_moment_z / (density * area));
};

using elements::tensor::as_vector;

int main(int argc, char **argv) {

  // MPI/PETSc/cpppetsc must be initialized before using it.
  const auto context = Context(&argc, &argv);

  // The beam with defined length is discretized into 10 segments
  const auto geometry = elements::mesh::refine_segment_mesh(
      elements::mesh::Mesh<Point>{
          {{0, 1}},
          {{0., 0.}, {beam_length, 0.}},
      },
      beam_length / 10);

  // Let's create the mesh and an assembler.
  const auto mesh = Mesh::fromConnectivity(
      coordinate_dimension, geometry.connectivity(),
      geometry.number_of_positions(), Element::degrees_of_freedom(), 0);

  auto assembler = Assembler();

  Properties properties = {
      young_modulus, shear_modulus, shear_correction_factor_y,
      area,          area_moment_z,
  };

  for (const auto &element : mesh.localElements()) {
    elements::tensor::Tensor<Mesh::real_type, coordinate_dimension>
        element_axis;
    as_vector(&element_axis) =
        as_vector(&geometry.position_of_vertex(
            geometry.connectivity().at(element.index()).at(1))) -
        as_vector(&geometry.position_of_vertex(
            geometry.connectivity().at(element.index()).at(0)));

    assembler.emplaceElement(element,
                             Element::Element(timoshenko_beam_stiffness_matrix(
                                 element_axis, properties)),
                             timoshenko_beam_consistent_mass_matrix(
                                 element_axis, properties, density));
  }

  // Let's assemble the stiffness and mass matrices.
  const auto K = [&]() {
    auto K = Mesh::matrix_type::fromMesh(mesh);
    assembler.assembleStiffnessMatrix(Mesh::vector_type::fromLocalMesh(mesh),
                                      Element::Time{0.}, &K);
    K.finalize();
    return K;
  }();

  const auto M = [&]() {
    auto M = Mesh::matrix_type::fromMesh(mesh);
    assembler.assembleMassMatrix(&M);
    M.finalize();
    return M;
  }();

  // Finally, we can compute the eigenvalues.
  const auto number_of_eigenvalues = Mesh::size_type(10);
  const auto result = cppslepc::computeSmallestEigenvalues(
      cpppetsc::computeElementsOfMatrix(K),
      cpppetsc::computeElementsOfMatrix(M), number_of_eigenvalues);

  // Let's compare the results to the analytical solution.
  for (const auto [n, value] : ranges::views::enumerate(result | rv::drop(2)) |
                                   rv::drop(1) | rv::take(4))
    Policy::handleError(
        PetscPrintf(Policy::communicator(),
                    "natural frequency %zu: %.4f (error: %3.2f%%)\n",
                    std::size_t{n}, std::sqrt(value).real(),
                    (std::sqrt(value).real() - analytic_result(n)) /
                        analytic_result(n) * 100));
}
