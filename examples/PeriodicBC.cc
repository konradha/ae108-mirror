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
#include "ae108/cpppetsc/Context.h"
#include "ae108/cpppetsc/GeneralizedMeshBoundaryCondition.h"
#include "ae108/cpppetsc/Mesh.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cpppetsc/createTransformInput.h"
#include "ae108/elements/CoreElement.h"
#include "ae108/elements/embedding/IsoparametricEmbedding.h"
#include "ae108/elements/integrator/IsoparametricIntegrator.h"
#include "ae108/elements/materialmodels/Hookean.h"
#include "ae108/elements/quadrature/Quadrature.h"
#include "ae108/elements/shape/Tri3.h"
#include "ae108/solve/TransformingSolver.h"
#include "ae108/solve/boundaryConditionsToTransform.h"
#include <array>

using namespace ae108;

using Policy = cpppetsc::ParallelComputePolicy;
using Context = cpppetsc::Context<Policy>;
using Mesh = cpppetsc::Mesh<Policy>;
using Vector = cpppetsc::Vector<Policy>;
using BoundaryCondition = cpppetsc::GeneralizedMeshBoundaryCondition<Mesh>;

namespace embedding = ae108::elements::embedding;
namespace integrator = ae108::elements::integrator;
namespace materialmodels = ae108::elements::materialmodels;
namespace quadrature = ae108::elements::quadrature;
namespace shape = ae108::elements::shape;

// In this example we will simulate two elements, A and B.
// Each of these elements has three vertices.

//  3 ------ 2
//  |  \  B  |
//  |  A  \  |
//  0--------1

// Let's specify the parameters of this mesh.

constexpr auto number_of_vertices_per_element = Mesh::size_type{3};
constexpr auto number_of_elements = Mesh::size_type{2};
constexpr auto number_of_vertices = Mesh::size_type{4};
constexpr auto dimension = Mesh::size_type{2};
constexpr auto dof_per_vertex = Mesh::size_type{2};

// The connectivity specifies the vertex indices per Element.

using Connectivity =
    std::array<std::array<Mesh::size_type, number_of_vertices_per_element>,
               number_of_elements>;
constexpr auto connectivity = Connectivity{{
    {{0, 1, 3}}, // vertices of element A
    {{1, 2, 3}}, // vertices of element B
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
}};

// Let's specify how the elements behave when deformed.

// First, we select a Hookean (linear elastic) material model in 2D.

using MaterialModel =
    materialmodels::Hookean<dimension, Vector::value_type, Vector::real_type>;

// We'll choose the Tri3 shape functions (3 shape functions) that are defined
// in 2D reference space.

using Shape = shape::Tri3;

// Since also our physical space is 2D we may use the isoparametric embedding.

using Embedding = embedding::IsoparametricEmbedding<Shape>;

// Let's choose a simple quadrature rule in reference space.

using Quadrature = quadrature::Quadrature<quadrature::QuadratureType::Simplex,
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

using Solver = solve::TransformingSolver<Assembler>;

int main(int argc, char **argv) {
  // MPI/PETSc/cpppetsc must be initialized before using it.

  const auto context = Context(&argc, &argv);

  // Let's create the mesh, an assembler, and a material model.

  const auto mesh = Mesh::fromConnectivity(dimension, connectivity,
                                           number_of_vertices, dof_per_vertex);
  auto assembler = Assembler();
  const auto model = MaterialModel(1.0, .2);

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
        }})));
  }

  // We need to create a solver. We do not use the time, so we can set it to
  // zero.

  const auto solver = Solver(&mesh);
  const auto time = Element::Time{0.};

  // Before we can produce meaningful results, we need to specify boundary
  // conditions. Let's fix the nodes at x=0 and pull on the nodes at x=2.

  // Note that we are using the degrees of freedom of vertex 3 to specify the
  // boundary conditions at vertex 2.

  std::vector<BoundaryCondition> boundary_conditions;
  for (const auto &vertex : mesh.localVertices()) {
    switch (vertex.index()) {
    case 0: {
      // The displacement in x direction is zero.
      boundary_conditions.push_back({{vertex.index(), 0}, {}, 0.});
      // The displacement in y direction is zero.
      boundary_conditions.push_back({{vertex.index(), 1}, {}, 0.});
      break;
    }
    case 1: {
      // The displacement in x direction is .5.
      boundary_conditions.push_back({{vertex.index(), 0}, {}, .5});
      // The displacement in y direction is zero.
      boundary_conditions.push_back({{vertex.index(), 1}, {}, 0.});
      break;
    }
    case 2: {
      const auto source_vertex = Mesh::size_type{3};
      // The displacement in x direction is .5 plus the displacement in
      // x direction of the vertex with index 3.
      boundary_conditions.push_back({{vertex.index(), 0},
                                     {
                                         {1., {source_vertex, 0}},
                                     },
                                     .5});
      // The displacement in y direction is the displacement in y direction of
      // the vertex with index 3.
      boundary_conditions.push_back({{vertex.index(), 1},
                                     {
                                         {1., {source_vertex, 1}},
                                     },
                                     0.});
      break;
    }
    }
  }

  // We are ready to minimize the energy.

  const auto transform =
      solve::boundaryConditionsToTransform(boundary_conditions, mesh);

  const auto result = apply(
      transform,
      solver.computeSolution(transform, createTransformInput(transform.matrix),
                             time, &assembler));

  // As a final step, we print the result after reordering to "canonical"
  // order (i.e. [ vertex-0-dof-0, vertex-0-dof-1, vertex-1-dof-0,
  // vertex-1-dof-1, ...]).

  const auto global_result =
      Vector::fromDistributedInCanonicalOrder(result, mesh);
  if (Policy::isPrimaryRank())
    global_result.unwrap().print();
}
