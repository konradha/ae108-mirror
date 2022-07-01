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

#include "ae108/assembly/Assembler.h"
#include "ae108/cpppetsc/Context.h"
#include "ae108/cpppetsc/Mesh.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cpppetsc/Viewer.h"
#include "ae108/cpppetsc/createVectorFromSource.h"
#include "ae108/cpppetsc/setName.h"
#include "ae108/cpppetsc/writeToViewer.h"
#include "ae108/elements/TimoshenkoBeamElement.h"
#include "ae108/elements/tensor/as_vector.h"
#include "ae108/solve/NonlinearSolver.h"

using namespace ae108;

using Policy = cpppetsc::ParallelComputePolicy;
using Context = cpppetsc::Context<Policy>;
using Mesh = cpppetsc::Mesh<Policy>;
using Vector = cpppetsc::Vector<Policy>;
using Viewer = cpppetsc::Viewer<Policy>;
using BoundaryCondition = cpppetsc::MeshBoundaryCondition<Mesh>;

// In this example we will simulate six elements, A, B, C, D, E, F.
// Each of these elements has two vertices.

//  4 --D---3
//    \     |  |
//      F   E    C
//        \ |      |
//  0---A---1---B---2

// Let's specify the parameters of this mesh.

constexpr auto number_of_vertices_per_element = Mesh::size_type{2};
constexpr auto number_of_elements = Mesh::size_type{6};
constexpr auto number_of_vertices = Mesh::size_type{5};
constexpr auto coordinate_dimension = Mesh::size_type{2};
constexpr auto dof_per_vertex = Mesh::size_type{3};
constexpr auto dof_per_element = Mesh::size_type{0};

// The connectivity specifies the vertex indices per Element.

using Connectivity =
    std::array<std::array<Mesh::size_type, number_of_vertices_per_element>,
               number_of_elements>;
constexpr auto connectivity = Connectivity{{
    {{0, 1}}, // vertices of element A
    {{1, 2}}, // vertices of element B
    {{2, 3}}, // vertices of element C
    {{3, 4}}, // vertices of element D
    {{1, 3}}, // vertices of element E
    {{1, 4}}, // vertices of element F
}};

// Each of the vertices is located at the following coordinates in physical
// space.

using VertexPositions =
    std::array<std::array<Mesh::real_type, coordinate_dimension>,
               number_of_vertices>;
constexpr auto vertex_positions = VertexPositions{{
    {{0., 0.}},
    {{1., 0.}},
    {{2., 0.}},
    {{1., 1.}},
    {{0., 1.}},
}};

// Now we are ready to select the Timoshenko beam element

using Element =
    elements::TimoshenkoBeamElement<coordinate_dimension, Vector::value_type,
                                    Vector::real_type>;

// The Timoshenko beam element comes with a number of geometrical and material
// related properties. These are stored in the Properties struct

using Properties =
    elements::TimoshenkoBeamProperties<Mesh::real_type, coordinate_dimension>;

// Let us define some parameters the linear elastic beam of rectangular cross
// section

constexpr Mesh::real_type young_modulus = 1.;
constexpr Mesh::real_type poisson_ratio = 0.3;
constexpr Mesh::real_type shear_modulus =
    young_modulus / (2 * (1 + poisson_ratio));
constexpr Mesh::real_type shear_correction_factor_y = 1.2;
constexpr Mesh::real_type thickness = 1.;
constexpr Mesh::real_type width = 0.1;
constexpr Mesh::real_type area = width * thickness;
constexpr Mesh::real_type area_moment_z = thickness * std::pow(width, 3) / 12.;

// We will assemble e.g. energy using a collection of elements. This is done by
// the assembler. (The list DefaultFeaturePlugins contain the features (e.g.
// energy) that most elements support.)

using Assembler =
    assembly::Assembler<Element, assembly::DefaultFeaturePlugins, Policy>;

// Our goals is to minimize the energy. This is done by the solver.

using Solver = solve::NonlinearSolver<Assembler>;

int main(int argc, char **argv) {
  // MPI/PETSc/cpppetsc must be initialized before using it.

  const auto context = Context(&argc, &argv);

  // Let's create the mesh and an assembler.

  const auto mesh =
      Mesh::fromConnectivity(coordinate_dimension, connectivity,
                             number_of_vertices, dof_per_vertex, 0);
  auto assembler = Assembler();

  Properties properties = {young_modulus, shear_modulus,
                           shear_correction_factor_y, area, area_moment_z};

  // Depending on whether we use MPI, our mesh may be distributed and not all
  // elements are present on this computational node.

  // Let's add those elements that are "local" to the assembler.

  for (const auto &element : mesh.localElements()) {
    elements::tensor::Tensor<Mesh::real_type, coordinate_dimension>
        element_axis;
    elements::tensor::as_vector(&element_axis) =
        elements::tensor::as_vector(
            &vertex_positions.at(connectivity.at(element.index()).at(1))) -
        elements::tensor::as_vector(
            &vertex_positions.at(connectivity.at(element.index()).at(0)));

    assembler.emplaceElement(
        element, timoshenko_beam_stiffness_matrix(element_axis, properties));
  }

  // We need to create a solver. We do not use the time, so we can set it to
  // zero.

  const auto solver = Solver(&mesh);
  const auto time = Element::Time{0.};

  // Before we can produce meaningful results, we need to specify boundary
  // conditions. Let's fix the nodes at x=0 and pull on the nodes at x=2.

  std::vector<BoundaryCondition> boundary_conditions;
  for (const auto &vertex : mesh.localVertices()) {
    switch (vertex.index()) {
    case 0:
    case 4: {
      // The displacement in x direction is zero.
      boundary_conditions.push_back({vertex, 0, 0.});
      // The displacement in y direction is zero.
      boundary_conditions.push_back({vertex, 1, 0.});
      break;
    }
    case 2: {
      // The displacement in y direction is -0.1
      boundary_conditions.push_back({vertex, 1, -0.1});
      break;
    }
    }
  }

  // We are ready to minimize the energy.

  auto result = solver.computeSolution(
      boundary_conditions, Vector::fromGlobalMesh(mesh), time, &assembler);

  // Let us now print the result.
  // We could use: result.unwrap().print();

  // However, the output of this command changes when running the application
  // in MPI-parallel mode. This happens because the vector "result" is
  // distributed between the ranks in this case.

  // So we gather all the results locally in the canonical data format first:
  // [ vertex-0-dof-0, vertex-0-dof-1, vertex-1-dof-0, vertex-1-dof-1, ...].

  const auto global_result =
      Vector::fromDistributedInCanonicalOrder(result, mesh);

  if (Policy::isPrimaryRank())
    global_result.unwrap().print();

  // As a final step, we write the result to a file.

  // We may view this result in e.g. ParaView by generating an XDMF file with
  // the `generate_xdmf.py` script in the repository, and opening this file
  // with ParaView ("XDMF Reader").

  // First we collect the coordinates in a vector.
  using DataSource = std::function<void(Mesh::size_type, Mesh::value_type *)>;
  auto coordinates = cpppetsc::createVectorFromSource(
      mesh, coordinate_dimension,
      DataSource([&](const Mesh::size_type index, Mesh::value_type *const out) {
        const auto &position = vertex_positions.at(index);
        std::copy(position.begin(), position.end(), out);
      }));
  cpppetsc::setName("coordinates", &coordinates);

  // Now we write the mesh to a file.
  auto viewer =
      Viewer::fromHdf5FilePath("timoshenko_beam.ae108", Viewer::Mode::write);
  cpppetsc::writeToViewer(mesh, coordinates, &viewer);

  // Let's write the result to the file.
  cpppetsc::setName("result", &result);
  cpppetsc::writeToViewer(result, &viewer);

  fprintf(stderr, "The data has been written to the file "
                  "\"timoshenko_beam.ae108\".\n");
}
