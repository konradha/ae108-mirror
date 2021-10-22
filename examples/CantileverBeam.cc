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

#include "ae108/assembly/Assembler.h"
#include "ae108/assembly/AssemblerGroup.h"
#include "ae108/cpppetsc/Context.h"
#include "ae108/cpppetsc/Mesh.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/elements/ForceElement.h"
#include "ae108/elements/TimoshenkoBeamElement.h"
#include "ae108/elements/tensor/as_vector.h"
#include "ae108/solve/NonlinearSolver.h"

using namespace ae108;

using Policy = cpppetsc::ParallelComputePolicy;
using Context = cpppetsc::Context<Policy>;
using Mesh = cpppetsc::Mesh<Policy>;
using Vector = cpppetsc::Vector<Policy>;
using BoundaryCondition = cpppetsc::MeshBoundaryCondition<Mesh>;

// In this example we will calculate the tip displacement of a cantilver beam
// 'A' with a force 'B' on the free end.

//        B ^
//          |
//  0---A---1

// Let's specify the parameters of this mesh.

constexpr auto number_of_vertices_per_element = Mesh::size_type{2};
constexpr auto number_of_elements = Mesh::size_type{2};
constexpr auto number_of_vertices = Mesh::size_type{2};
constexpr auto coordinate_dimension = Mesh::CoordinateDimension{3};
constexpr auto topological_dimension = Mesh::TopologicalDimension{1};
constexpr auto dof_per_vertex = Mesh::size_type{6};
constexpr auto dof_per_element = Mesh::size_type{0};

// The connectivity specifies the vertex indices per Element.

using Connectivity =
    std::array<std::array<Mesh::size_type, number_of_vertices_per_element>,
               number_of_elements>;
constexpr auto connectivity = Connectivity{{
    {{0, 1}}, // vertices of element A
    {{1}},    // vertex of tip load B
}};

// Vertices 1 and 2 are located at the following coordinates in physical
// space.

using VertexPositions =
    std::array<std::array<Mesh::real_type, coordinate_dimension>,
               number_of_vertices>;
constexpr auto vertex_positions = VertexPositions{{
    {{0., 0., 0.}},
    {{1., 0., 0.}},
}};

// Now we are ready to select the Timoshenko beam element

using Element =
    elements::TimoshenkoBeamElement<coordinate_dimension, Vector::value_type,
                                    Vector::real_type>;

// The Timoshenko beam element comes with a number of geometrical and material
// related properties. These are stored in the Properties struct

using Properties =
    elements::TimoshenkoBeamProperties<Mesh::real_type, coordinate_dimension>;

// Let us define some parameters for the linear elastic beam of circular cross
// section

constexpr Mesh::real_type young_modulus = 1.;
constexpr Mesh::real_type poisson_ratio = 0.3;
constexpr Mesh::real_type shear_modulus =
    young_modulus / (2 * (1 + poisson_ratio));
constexpr Mesh::real_type shear_correction_factor =
    (7 + 6 * poisson_ratio) / 6 / (1 + poisson_ratio); // Cowper (1966)
constexpr Mesh::real_type radius = 0.05;
constexpr Mesh::real_type area = radius * radius * M_PI;
constexpr Mesh::real_type area_moment =
    M_PI_4 * radius * radius * radius * radius;
constexpr Mesh::real_type polar_moment =
    M_PI_2 * radius * radius * radius * radius;

// We will assemble e.g. energy using a collection of elements. This is done by
// the assembler. (The list DefaultFeaturePlugins contain the features (e.g.
// energy) that most elements support.)

using Assembler =
    assembly::Assembler<Element, assembly::DefaultFeaturePlugins, Policy>;

using ForceElement = elements::ForceElement<dof_per_vertex, Vector::value_type,
                                            Vector::real_type>;
using ForceAssembler =
    assembly::Assembler<ForceElement, assembly::DefaultFeaturePlugins, Policy>;

// Finally we combine the assemblers in a AssemblerGroup.

using GroupAssembler = assembly::AssemblerGroup<Assembler, ForceAssembler>;

// Our goals is to minimize the energy. This is done by the solver.

using Solver = solve::NonlinearSolver<GroupAssembler>;

int main(int argc, char **argv) {
  // MPI/PETSc/cpppetsc must be initialized before using it.

  const auto context = Context(&argc, &argv);

  // Let's create the mesh and an assembler.

  const auto mesh = Mesh::fromConnectivity(
      topological_dimension, coordinate_dimension, connectivity,
      number_of_vertices, dof_per_vertex, 0);
  auto assembler = GroupAssembler();
  auto &element_assembler = assembler.get<0>();
  auto &force_assembler = assembler.get<1>();

  Properties properties = {young_modulus,
                           shear_modulus,
                           shear_correction_factor,
                           shear_correction_factor,
                           area,
                           area_moment,
                           area_moment,
                           polar_moment};

  // Depending on whether we use MPI, our mesh may be distributed and not all
  // elements are present on this computational node.

  // Let's add those elements that are "local" to the assembler.

  for (const auto &element : mesh.localElements()) {
    switch (element.index()) {
    case 0: {
      elements::tensor::Tensor<Mesh::real_type, coordinate_dimension>
          element_axis;
      elements::tensor::as_vector(&element_axis) =
          elements::tensor::as_vector(
              &vertex_positions.at(connectivity.at(element.index()).at(1))) -
          elements::tensor::as_vector(
              &vertex_positions.at(connectivity.at(element.index()).at(0)));

      element_assembler.emplaceElement(
          element, timoshenko_beam_stiffness_matrix(element_axis, properties));
      break;
    }
    case 1: {
      force_assembler.emplaceElement(
          element, ForceElement::Force{{0., 1e-6, 0., 0., 0., 0.}});
      break;
    }
    }
  }
  // We need to create a solver. We do not use the time, so we can set it to
  // zero.

  const auto solver = Solver(&mesh);
  const auto time = Element::Time{0.};

  // Before we can produce meaningful results, we need to specify boundary
  // conditions. Let's clamp the beam completely (i.e. all dofs) at (0,0)

  std::vector<BoundaryCondition> boundary_conditions;
  for (const auto &vertex : mesh.localVertices()) {
    switch (vertex.index()) {
    case 0: {
      // The displacement in x direction is zero.
      boundary_conditions.push_back({vertex, 0, 0.});
      // The displacement in y direction is zero.
      boundary_conditions.push_back({vertex, 1, 0.});
      // The displacement in z direction is zero.
      boundary_conditions.push_back({vertex, 2, 0.});
      // The rotation around x axis is zero.
      boundary_conditions.push_back({vertex, 3, 0.});
      // The rotation around y axis is zero.
      boundary_conditions.push_back({vertex, 4, 0.});
      // The rotation around z axis is zero.
      boundary_conditions.push_back({vertex, 5, 0.});
      break;
    }
    }
  }

  // We are ready to minimize the energy.

  auto result = solver.computeSolution(
      boundary_conditions, Vector::fromGlobalMesh(mesh), time, &assembler);

  // Let us now print the result.
  // We could use: result.unwrap().print();

  // However, the output of this command changes when running the
  // application in MPI-parallel mode. This happens because the vector
  // "result" is distributed between the ranks in this case.

  // So we gather all the results locally in the canonical data format
  // first: [ vertex-0-dof-0, vertex-0-dof-1, vertex-1-dof-0,
  // vertex-1-dof-1, ...].

  const auto global_result =
      Vector::fromDistributedInCanonicalOrder(result, mesh);

  if (Policy::isPrimaryRank())
    global_result.unwrap().print();

  PetscPrin
}
