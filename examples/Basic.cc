// © 2020 ETH Zurich, Mechanics and Materials Lab
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

#include "ae108/assembly/Assembler.h"
#include "ae108/cpppetsc/Context.h"
#include "ae108/cpppetsc/Mesh.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/elements/CoreElement.h"
#include "ae108/elements/embedding/IsoparametricEmbedding.h"
#include "ae108/elements/integrator/IsoparametricIntegrator.h"
#include "ae108/elements/materialmodels/Hookean.h"
#include "ae108/elements/quadrature/Quadrature.h"
#include "ae108/elements/shape/Quad4.h"
#include "ae108/solve/NonlinearSolver.h"
#include <array>

using namespace ae108;

using Policy = cpppetsc::ParallelComputePolicy;
using Context = cpppetsc::Context<Policy>;
using Mesh = cpppetsc::Mesh<Policy>;
using Vector = cpppetsc::Vector<Policy>;
using BoundaryCondition = cpppetsc::MeshBoundaryCondition<Mesh>;

namespace embedding = ae108::elements::embedding;
namespace integrator = ae108::elements::integrator;
namespace materialmodels = ae108::elements::materialmodels;
namespace quadrature = ae108::elements::quadrature;
namespace shape = ae108::elements::shape;

// In this example we will simulate two elements, A and B.
// Each of these elements has four vertices.

//  3 ------- 2 ------- 5
//  |         |         |
//  |    A    |    B    |
//  0---------1---------4

// Let's specify the parameters of this mesh.

constexpr auto number_of_vertices_per_element = Mesh::size_type{4};
constexpr auto number_of_elements = Mesh::size_type{2};
constexpr auto number_of_vertices = Mesh::size_type{6};
constexpr auto dimension = Mesh::size_type{2};
constexpr auto dof_per_vertex = Mesh::size_type{2};

// The connectivity specifies the vertex indices per Element.

using Connectivity =
    std::array<std::array<Mesh::size_type, number_of_vertices_per_element>,
               number_of_elements>;
constexpr auto connectivity = Connectivity{{
    {{0, 1, 2, 3}}, // vertices of element A
    {{1, 4, 5, 2}}, // vertices of element B
}};

// Each of the vertices is located at the following coordinates in physical
// space.

using VertexPositions =
    std::array<std::array<Mesh::real_type, dimension>, number_of_vertices>;
constexpr auto vertex_positions = VertexPositions{{
    {{0., 0.}},
    {{1., 0.}},
    {{1., 1.}},
    {{0., 1.}},
    {{2., 0.}},
    {{2., 1.}},
}};

// Let's specify how the elements behave when deformed.

// First, we select a Hookean (linear elastic) material model in 2D.

using MaterialModel =
    materialmodels::Hookean<dimension, Vector::value_type, Vector::real_type>;

// We'll choose the Quad4 shape functions (4 shape functions) that are defined
// in 2D reference space.

using Shape = shape::Quad4;

// Since also our physical space is 2D we may use the isoparametric embedding.

using Embedding = embedding::IsoparametricEmbedding<Shape>;

// Let's choose a simple quadrature rule in reference space.

using Quadrature = quadrature::Quadrature<quadrature::QuadratureType::Cube,
                                          dimension, 1 /* order */>;

// We use an "Integrator" to integrate in physical space. Of course, it depends
// on the selected quadrature rule and the embedding. In the case of the
// isoparametric embedding, this embedding is specified by the "Shape".

using Integrator =
    integrator::IsoparametricIntegrator<Shape, Quadrature, Vector::value_type,
                                        Vector::real_type>;

// Now we are ready to select our element type.

using Element = elements::CoreElement<MaterialModel, Integrator,
                                      Vector::value_type, Vector::real_type>;

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

  // Let's create the mesh, an assembler, and a material model.

  const auto mesh = Mesh::fromConnectivity(dimension, connectivity,
                                           number_of_vertices, dof_per_vertex);
  auto assembler = Assembler();
  const auto model = MaterialModel(1.0, 0.);

  // Depending on whether we use MPI, our mesh may be distributed and not all
  // elements are present on this computational node.

  // Let's add those elements that are "local" to the assembler.

  for (const auto &element : mesh.localElements()) {
    assembler.emplaceElement(
        element, model,
        Integrator(Embedding(Embedding::Collection<Embedding::PhysicalPoint>{{
            vertex_positions.at(connectivity.at(element.index()).at(0)),
            vertex_positions.at(connectivity.at(element.index()).at(1)),
            vertex_positions.at(connectivity.at(element.index()).at(2)),
            vertex_positions.at(connectivity.at(element.index()).at(3)),
        }})));
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
    case 3: {
      // The displacement in x direction is zero.
      boundary_conditions.push_back({vertex, 0, 0.});
      // The displacement in y direction is zero.
      boundary_conditions.push_back({vertex, 1, 0.});
      break;
    }
    case 4:
    case 5: {
      // The displacement in x direction is .5.
      boundary_conditions.push_back({vertex, 0, .5});
      // The displacement in y direction is 0.
      boundary_conditions.push_back({vertex, 1, 0.});
      break;
    }
    }
  }

  // We are ready to minimize the energy.

  const auto result = solver.computeSolution(
      boundary_conditions, Vector::fromGlobalMesh(mesh), time, &assembler);

  // As a final step, we print the result.
  // We could use: result.unwrap().print();

  // However, the output of this command changes when running the application
  // in MPI-parallel mode. This happens because the vector "result" is
  // distributed between the ranks in this case.

  // So we gather all the results locally in the canonical data format first:
  // [ vertex-0-dof-0, vertex-0-dof-1, vertex-1-dof-0, vertex-1-dof-1, ...].

  const auto global_result =
      Vector::fromDistributedInCanonicalOrder(result, mesh);

  // Then we print this global vector only on the primary rank.
  // We expect the degrees of freedom in x direction to be between 0 and .5,
  // and the degrees of freedom in y direction to be 0:
  // Vertex 0: (0.0, 0.0)
  // Vertex 1: (.25, 0.0)
  // Vertex 2: (.25, 0.0)
  // Vertex 3: (0.0, 0.0)
  // Vertex 4: (0.5, 0.0)
  // Vertex 5: (0.5, 0.0)

  if (Policy::isPrimaryRank())
    global_result.unwrap().print();
}
