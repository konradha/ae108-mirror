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
#include "ae108/cpppetsc/GeneralizedMeshBoundaryCondition.h"
#include "ae108/cpppetsc/Mesh.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cpppetsc/asThAT.h"
#include "ae108/cpppetsc/computeElementsOfMatrix.h"
#include "ae108/cppslepc/Context.h"
#include "ae108/cppslepc/computeSmallestEigenvalues.h"
#include "ae108/elements/ElementWithMass.h"
#include "ae108/elements/TimoshenkoBeamElement.h"
#include "ae108/elements/compute_mass_matrix.h"
#include "ae108/elements/mesh/refine_segment_mesh.h"
#include "ae108/elements/tensor/as_vector.h"
#include "ae108/solve/boundaryConditionsToTransform.h"
#include <array>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/take.hpp>

using namespace ae108;
namespace rv = ranges::cpp20::views;

using Policy = cpppetsc::ParallelComputePolicy;
using Context = cppslepc::Context<Policy>;
using Mesh = cpppetsc::Mesh<Policy>;
using Vector = cpppetsc::Vector<Policy>;
using BoundaryCondition = cpppetsc::GeneralizedMeshBoundaryCondition<Mesh>;

// In this example we investigate the plane wave propagation in an infinite
// two-dimensional square lattice as presented in "Wave propagation in
// two-dimensional periodic lattices", Phani et al. (2006), Section V.D
// https://doi.org/10.1121/1.2179748

//  3 --C---2
//  |       |
//  D       B
//  |       |
//  0---A---1

// Let's specify the parameters of this problem.
constexpr auto coordinate_dimension = Mesh::size_type{2};

using Point = std::array<Mesh::real_type, coordinate_dimension>;

// The element needs to have a ComputeMassMatrixTrait
using TimoshenkoBeamElement =
    elements::TimoshenkoBeamElement<coordinate_dimension, Vector::value_type,
                                    Vector::real_type>;
using Element = elements::ElementWithMass<TimoshenkoBeamElement>;

// The Timoshenko beam element comes with a number of geometrical and material
// related properties. These are stored in the Properties struct
using Properties =
    elements::TimoshenkoBeamProperties<Mesh::real_type, coordinate_dimension>;

// Make sure to add the AssembleMassMatrixPlugin, as we need the mass matrix in
// the eigenvalue problem
using Plugins =
    assembly::FeaturePlugins<assembly::plugins::AssembleEnergyPlugin,
                             assembly::plugins::AssembleForceVectorPlugin,
                             assembly::plugins::AssembleStiffnessMatrixPlugin,
                             assembly::plugins::AssembleMassMatrixPlugin>;

using Assembler = assembly::Assembler<Element, Plugins, Policy>;

// Let us define some parameters of the linear elastic beam of rectangular cross
// section as presented by Phani et al. (2009)
constexpr Mesh::real_type young_modulus = 1.;
constexpr Mesh::real_type poisson_ratio = 0.3;
constexpr Mesh::real_type shear_modulus =
    young_modulus / (2 * (1 + poisson_ratio));
constexpr Mesh::real_type shear_correction_factor_y = 1.2;
constexpr Mesh::real_type thickness = 1.;
constexpr Mesh::real_type slenderness = 50.;
constexpr Mesh::real_type beam_length = 1.;
constexpr Mesh::real_type width = 2. * 1.732050808 * beam_length / slenderness;
constexpr Mesh::real_type area = width * thickness;
constexpr Mesh::real_type density = 1.;
constexpr Mesh::real_type area_moment_z =
    thickness * width * width * width / 12.;

// resonance frequencies for a pinned-pinned beam
const auto analytic_result = [](std::size_t n) {
  return std::pow(n * M_PI / beam_length, 2) *
         std::sqrt(young_modulus * area_moment_z / (density * area));
};

using elements::tensor::as_vector;

int main(int argc, char **argv) {

  // MPI/PETSc/SLEPc must be initialized before using it.
  const auto context = Context(&argc, &argv);

  // Introduce the example.
  Policy::handleError(
      PetscPrintf(Policy::communicator(), "%s", std::string(90, '#').c_str()));
  Policy::handleError(PetscPrintf(
      Policy::communicator(),
      "\nCompare to Phani et al. (2006), https://doi.org/10.1121/1.2179748"
      "\nSection V.D., slenderness 50, Figure 14\n"));
  Policy::handleError(
      PetscPrintf(Policy::communicator(), "%s", std::string(90, '#').c_str()));
  Policy::handleError(PetscPrintf(Policy::communicator(), "\n"));
  Policy::handleError(
      PetscPrintf(Policy::communicator(), "%s", std::string(12, ' ').c_str()));

  // We define the square lattice (Phani et al., case D).
  const elements::mesh::Mesh<Point> geometry{
      {{0, 1}, {1, 2}, {2, 3}, {3, 0}},
      {{0., 0.}, {1., 0.}, {1., 1.}, {0., 1.}},
  };

  // We choose to discretize each side into 20 elements.
  const auto number_of_elements_per_beam = 20;
  const auto refined_geometry = elements::mesh::refine_segment_mesh(
      geometry, beam_length / number_of_elements_per_beam);

  // Let's create the mesh and an assembler.
  const auto mesh = Mesh::fromConnectivity(
      coordinate_dimension, refined_geometry.connectivity(),
      refined_geometry.number_of_positions(), Element::degrees_of_freedom(), 0);

  auto assembler = Assembler();

  const auto properties = Properties{
      young_modulus, shear_modulus, shear_correction_factor_y,
      area,          area_moment_z,
  };

  // Let's add those elements that are "local" to the assembler.
  for (const auto &element : mesh.localElements()) {
    elements::tensor::Tensor<Mesh::real_type, coordinate_dimension>
        element_axis;
    as_vector(&element_axis) =
        as_vector(&refined_geometry.position_of_vertex(
            refined_geometry.connectivity().at(element.index()).at(1))) -
        as_vector(&refined_geometry.position_of_vertex(
            refined_geometry.connectivity().at(element.index()).at(0)));

    // There is a constant factor of 0.5, because all beam segments lie on the
    // boundary of the unit cell and are hence shared with the respective
    // adjacent cell.
    assembler.emplaceElement(
        element,
        Element::Element(
            0.5 * timoshenko_beam_stiffness_matrix(element_axis, properties)),
        0.5 * timoshenko_beam_consistent_mass_matrix(element_axis, properties,
                                                     density));
  }

  // Construct stiffness and mass matrix.
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

  // Construct the periodicity relation in the vector source_of[target] (eq.
  // 11). If source_of[n] = m, node m depends on node n. If source_of[n] = n,
  // node n is independent.

  // Initialize all nodes as independent
  auto source_of =
      rv::iota(0, mesh.totalNumberOfVertices()) | ranges::to<std::vector>();
  // All corner nodes (1, 2, 3) are linked to the bottom left corner node (0).
  for (const auto i : rv::iota(0, 4))
    source_of[i] = 0;
  for (const auto i : rv::iota(0, number_of_elements_per_beam - 1)) {
    // All nodes on the top side C are linked to the bottom side A.
    source_of[4 + 2 * (number_of_elements_per_beam - 1) + i] =
        3 + 1 * (number_of_elements_per_beam - 1) - i;
    // All nodes on the right side B are linked to the left side D.
    source_of[4 + 3 * (number_of_elements_per_beam - 1) + i] =
        3 + 2 * (number_of_elements_per_beam - 1) - i;
  }

  for (const auto n : rv::iota(0, 8))
    Policy::handleError(
        PetscPrintf(Policy::communicator(), "Mode %d    ", int{n + 1}));

  // Define the high symmetry points of interest (see Phani et al.).
  const std::array<std::pair<char, Point>, 3> high_symmetry_points{{
      {'O', {0, 0}},
      {'A', {M_PI, 0}},
      {'B', {M_PI, M_PI}},
  }};

  // Solve for the eigenfrequencies at each k-point.
  for (const auto &point : high_symmetry_points) {
    // Construct the transformation matrix (eq. 12).
    const auto T = [&](const Point &wave_vector) {
      std::vector<BoundaryCondition> boundary_conditions;
      for (const auto &vertex : mesh.localVertices()) {
        const auto target = vertex.index();
        const auto source = source_of[target];
        if (source == target)
          continue;

        const auto k =
            PETSC_i *
            as_vector(&wave_vector)
                .dot(as_vector(&refined_geometry.position_of_vertex(target)) -
                     as_vector(&refined_geometry.position_of_vertex(source)));

        for (const auto dof : rv::iota(0, vertex.numberOfDofs())) {
          boundary_conditions.push_back({
              {target, dof},
              {{std::exp(k), {source, dof}}},
              0.,
          });
        }
      }
      return solve::boundaryConditionsToTransform(boundary_conditions, mesh)
          .matrix;
    }(point.second);

    // Transform mass and stiffness matrices (eq. 13).
    const auto KT = cpppetsc::asThAT(&K, &T);
    const auto MT = cpppetsc::asThAT(&M, &T);

    // Solve the generalized eigenvalue problem (eq. 14).
    const auto number_of_eigenvalues = Mesh::size_type(20);
    const auto result = cppslepc::computeSmallestEigenvalues(
        cpppetsc::computeElementsOfMatrix(KT),
        cpppetsc::computeElementsOfMatrix(MT), number_of_eigenvalues);

    // Present results.
    Policy::handleError(
        PetscPrintf(Policy::communicator(), "\nPoint %c", point.first));
    for (const auto value : result | rv::take(8))
      Policy::handleError(
          PetscPrintf(Policy::communicator(), "%10.2f",
                      std::sqrt(value).real() / analytic_result(1)));
  }
  Policy::handleError(PetscPrintf(Policy::communicator(), "\n"));
}