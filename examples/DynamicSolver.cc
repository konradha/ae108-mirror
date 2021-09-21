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

#include "ae108/solve/DynamicSolver.h"
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
#include "ae108/elements/shape/Seg2.h"
#include "ae108/solve/NonlinearSolver.h"
#include <array>
#include <iostream>

using namespace ae108;

using Policy = cpppetsc::ParallelComputePolicy;
using Context = cpppetsc::Context<Policy>;
using Mesh = cpppetsc::Mesh<Policy>;
using Vector = cpppetsc::Vector<Policy>;
using Matrix = cpppetsc::Matrix<Policy>;
using BoundaryCondition = cpppetsc::MeshBoundaryCondition<Mesh>;

namespace embedding = ae108::elements::embedding;
namespace integrator = ae108::elements::integrator;
namespace materialmodels = ae108::elements::materialmodels;
namespace quadrature = ae108::elements::quadrature;
namespace shape = ae108::elements::shape;

// In this example we will simulate one element, A.
// This element has two vertices.

//       A
//  0---------1

// Let's specify the parameters of this mesh.

constexpr auto number_of_vertices_per_element = Mesh::size_type{2};
constexpr auto number_of_elements = Mesh::size_type{1};
constexpr auto number_of_vertices = Mesh::size_type{2};
constexpr auto dimension = Mesh::size_type{1};
constexpr auto dof_per_vertex = Mesh::size_type{1};

// The connectivity specifies the vertex indices per element.

using Connectivity =
    std::array<std::array<Mesh::size_type, number_of_vertices_per_element>,
               number_of_elements>;
constexpr auto connectivity = Connectivity{{
    {{0, 1}}, // vertices of element A
}};

// Each of the vertices is located at the following coordinates in physical
// space.

using VertexPositions =
    std::array<std::array<Mesh::real_type, dimension>, number_of_vertices>;
constexpr auto vertex_positions = VertexPositions{{
    {{0.}},
    {{1.}},
}};

// Let's specify how the elements behave when deformed.

// First, we select a Hookean (linear elastic) material model in 1D.

using MaterialModel =
    materialmodels::Hookean<dimension, Vector::value_type, Vector::real_type>;

// We'll choose the Seg2 shape functions (2 shape functions) that are defined
// in 1D reference space.

using Shape = shape::Seg2;

// Since also our physical space is 1D we may use the isoparametric embedding.

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
// the assembler. (The list DefaultFeaturePlugins contain the features, e.g.
// energy, that most elements support.)

using Assembler =
    assembly::Assembler<Element, assembly::DefaultFeaturePlugins, Policy>;

// Our goal is to minimize the energy. This is done by the solver.

using Solver = solve::NonlinearSolver<Assembler>;

// In this example, we will also use a DynamicSolver to solve a dynamic problem.

using DynamicSolver = solve::DynamicSolver<Assembler>;

// For this purpose we specify Newmark parameters (case "average acceleration")
// and initial state.

const auto newmark_parameters = solve::dynamics::NewmarkParameters{.25, .5};
using State = DynamicSolver::state_type;

int main(int argc, char **argv) {
  // MPI/PETSc/cpppetsc must be initialized before using it.

  const auto context = Context(&argc, &argv);

  // Let's create the mesh, an assembler, and a material model.

  const auto mesh = Mesh::fromConnectivity(dimension, connectivity,
                                           number_of_vertices, dof_per_vertex);
  auto assembler = Assembler();
  const auto model = MaterialModel(1.0, 0.);

  // Depending on whether we use MPI, our mesh may be distributed and A
  // may not be present on this computational node.

  // Let's add those elements that are "local" to the assembler.

  for (const auto &element : mesh.localElements()) {
    assembler.emplaceElement(
        element, model,
        Integrator(Embedding(Embedding::Collection<Embedding::PhysicalPoint>{{
            vertex_positions.at(connectivity.at(element.index()).at(0)),
            vertex_positions.at(connectivity.at(element.index()).at(1)),
        }})));
  }

  // We need to create a solver.

  const auto solver = Solver(&mesh);

  // Before we can produce meaningful results, we need to specify boundary
  // conditions and initial state. Let's fix the nodes at x=0 and set initial
  // displacements of .5 at x=1.

  auto state = State{
      Vector::fromGlobalMesh(mesh), // displacements
      Vector::fromGlobalMesh(mesh), // velocities
      Vector::fromGlobalMesh(mesh), // accelerations
  };

  auto local_displacements = Vector::fromLocalMesh(mesh);
  std::vector<BoundaryCondition> boundary_conditions;
  for (const auto &vertex : mesh.localVertices()) {
    switch (vertex.index()) {
    case 0: {
      // The displacement in is zero.
      boundary_conditions.push_back({vertex, 0, 0.});
      break;
    }
    case 1: {
      // Initially, the displacement is .5
      vertex.setVertexData({.5}, &local_displacements);
      break;
    }
    }
  }
  mesh.copyToGlobalVector(local_displacements, &state.displacements);

  // Now we create a nonzero (lumped) mass matrix and a zero damping matrix.

  const auto mass = [&]() {
    auto mass = Matrix::fromMesh(mesh);
    auto replace = mass.preallocatedAssemblyView(1).replace();
    replace(0, 0) = .5;
    replace(1, 1) = .5;
    return mass;
  }();
  const auto damping = Matrix::fromMesh(mesh);

  // We are ready to run the simulation. We'll use a timestep of 1/10 und run
  // until t=10.

  const auto dynamic_solver = DynamicSolver{&mesh, &solver, newmark_parameters};
  const auto timestep = Element::Time{1. / 10.};

  for (auto time = Element::Time{0.}; time < Element::Time{10.};
       time += timestep) {
    state = dynamic_solver.computeSolution(boundary_conditions,
                                           std::move(state), time, timestep,
                                           mass, damping, &assembler);

    // Let's print the position of vertex 1 on the primary rank. We should see
    // it oscillate between extrema.

    const auto position =
        vertex_positions[1][0] +
        Vector::fromDistributedInCanonicalOrder(state.displacements, mesh)(1);

    if (Policy::isPrimaryRank())
      std::cout << std::real(position) << std::endl;
  }
}
