// ae108
// © 2023 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/elements/Bar.h"
#include "ae108/assembly/Assembler.h"
#include "ae108/cpppetsc/Context.h"
#include "ae108/cpppetsc/Mesh.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cpppetsc/Viewer.h"
#include "ae108/cpppetsc/createVectorFromSource.h"
#include "ae108/cpppetsc/setName.h"
#include "ae108/cpppetsc/writeToViewer.h"
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

//  3---C---2
//  | \   / |
//  D  E/F  B
//  | /   \ |
//  0---A---1

// Let's specify the parameters of this mesh.

constexpr auto number_of_vertices_per_element = Mesh::size_type{2};
constexpr auto number_of_elements = Mesh::size_type{6};
constexpr auto number_of_vertices = Mesh::size_type{4};
constexpr auto coordinate_dimension = Mesh::size_type{3};
constexpr auto dof_per_vertex = Mesh::size_type{coordinate_dimension};
constexpr auto dof_per_element = Mesh::size_type{0};

// The connectivity specifies the vertex indices per Element.

using Connectivity =
    std::array<std::array<Mesh::size_type, number_of_vertices_per_element>,
               number_of_elements>;
constexpr auto connectivity = Connectivity{{
    {{0, 1}}, // vertices of element A
    {{1, 2}}, // vertices of element B
    {{2, 3}}, // vertices of element C
    {{3, 0}}, // vertices of element D
    {{0, 2}}, // vertices of element E
    {{3, 1}}, // vertices of element F
}};

// Each of the vertices is located at the following coordinates in physical
// space.

using VertexPositions =
    std::array<std::array<Mesh::real_type, coordinate_dimension>,
               number_of_vertices>;
constexpr auto vertex_positions = VertexPositions{{
    {{0., 0., 0.}},
    {{1., 0., 0.}},
    {{1., 1., 0.}},
    {{0., 1., 0.}},
}};

// Now we are ready to select the Bar element.

using Element =
    elements::Bar<coordinate_dimension, Vector::value_type, Vector::real_type>;

// The bar element comes with geometrical and material
// properties. These are stored in the Properties struct.

using Properties = elements::BarProperties<Mesh::real_type>;

// Let us define some parameters for the bar.

constexpr Mesh::real_type young_modulus = 1.;
constexpr Mesh::real_type cross_section = 0.1;

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

  const auto properties = Properties{young_modulus, cross_section};

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
        element, elements::bar_stiffness_matrix<coordinate_dimension>(
                     element_axis, properties));
  }

  // We need to create a solver. We do not use the time, so we can set it to
  // zero.

  const auto solver = Solver(&mesh);
  const auto time = Element::Time{0.};

  // Before we can produce meaningful results, we need to specify boundary
  // conditions. Let's fix the node at (0, 0) and pull on the node at (1, 1).

  std::vector<BoundaryCondition> boundary_conditions;
  for (const auto &vertex : mesh.localVertices()) {
    // The displacement in z direction is zero.
    boundary_conditions.push_back({vertex, 2, 0.});
    switch (vertex.index()) {
    case 0: {
      // The displacement in x direction is zero.
      boundary_conditions.push_back({vertex, 0, 0.});
      // The displacement in y direction is zero.
      boundary_conditions.push_back({vertex, 1, 0.});
      break;
    }
    case 2: {
      // The displacement in x direction is 0.5.
      boundary_conditions.push_back({vertex, 0, 0.5});
      // The displacement in y direction is 0.5.
      boundary_conditions.push_back({vertex, 1, 0.5});
      break;
    }
    }
  }

  // We are ready to minimize the energy.

  auto result = solver.computeSolution(
      boundary_conditions, Vector::fromGlobalMesh(mesh), time, &assembler);

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

  auto viewer = Viewer::fromHdf5FilePath("bar.ae108", Viewer::Mode::write);
  cpppetsc::writeToViewer(mesh, coordinates, &viewer);

  // Let's write the result to the file.

  cpppetsc::setName("result", &result);
  cpppetsc::writeToViewer(result, &viewer);

  fprintf(stderr, "The data has been written to the file "
                  "\"bar.ae108\".\n");
}
