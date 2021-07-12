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

#include "ae108/elements/ForceElement.h"
#include "ae108/assembly/Assembler.h"
#include "ae108/assembly/AssemblerGroup.h"
#include "ae108/cpppetsc/Context.h"
#include "ae108/cpppetsc/Mesh.h"
#include "ae108/cpppetsc/ParallelComputePolicy.h"
#include "ae108/cpppetsc/Vector.h"
#include "ae108/cpppetsc/Viewer.h"
#include "ae108/cpppetsc/createVectorFromSource.h"
#include "ae108/cpppetsc/setName.h"
#include "ae108/cpppetsc/writeToViewer.h"
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

// In this example we will simulate an element A.
// A has four vertices.

//  3 ------- 2/B
//  |         |
//  |    A    |
//  0---------1

// In addition, we add forces at the x, y = 1. This is also modelled using a
// single vertex element at vertex 2.

// Let's specify the parameters of this mesh.

constexpr auto number_of_elements = Mesh::size_type{1 + 1};
constexpr auto number_of_vertices = Mesh::size_type{4};
constexpr auto dimension = Mesh::size_type{2};
constexpr auto dof_per_vertex = Mesh::size_type{2};

// The connectivity specifies the vertex indices per Element.

using Connectivity =
    std::array<std::vector<Mesh::size_type>, number_of_elements>;
const auto connectivity = Connectivity{{
    {{0, 1, 2, 3}}, // vertices of element A
    {{2}},          // vertex of force element B
}};

// Each of the vertices is located at the following coordinates in physical
// space.

using VertexPositions =
    std::array<std::array<Mesh::value_type, dimension>, number_of_vertices>;
constexpr auto vertex_positions = VertexPositions{{
    {{0., 0.}},
    {{1., 0.}},
    {{1., 1.}},
    {{0., 1.}},
}};

// Let's specify how the elements behave when deformed.

// First, we select a Hookean (linear elastic) material model in 2D.

using MaterialModel = materialmodels::Hookean<dimension>;

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

using Integrator = integrator::IsoparametricIntegrator<Shape, Quadrature>;

// Now we are ready to select our element type.

using Element = elements::CoreElement<MaterialModel, Integrator>;

// We will assemble e.g. energy using a collection of elements. This is done by
// the assembler. (The list DefaultFeaturePlugins contain the features (e.g.
// energy) that most elements support.)

using Assembler =
    assembly::Assembler<Element, assembly::DefaultFeaturePlugins, Policy>;

// In addition, we collect the force elements in another assembler.

using ForceElement = elements::ForceElement<dof_per_vertex>;
using ForceAssembler =
    assembly::Assembler<ForceElement, assembly::DefaultFeaturePlugins, Policy>;

// Finally we combine the assemblers in a AssemblerGroup.

using GroupAssembler = assembly::AssemblerGroup<Assembler, ForceAssembler>;

// Our goals is to minimize the energy. This is done by the solver.

using Solver = solve::NonlinearSolver<GroupAssembler>;

int main(int argc, char **argv) {
  // MPI/PETSc/cpppetsc must be initialized before using it.

  const auto context = Context(&argc, &argv);

  // Let's create the mesh, an assembler, and a material model.

  const auto mesh = Mesh::fromConnectivity(dimension, connectivity,
                                           number_of_vertices, dof_per_vertex);
  auto assembler = GroupAssembler();
  auto &element_assembler = assembler.get<0>();
  auto &force_assembler = assembler.get<1>();
  const auto model = MaterialModel(1.0, 0.0);

  // Depending on whether we use MPI, our mesh may be distributed and not all
  // elements are present on this computational node.

  // Let's add those elements that are "local" to the assembler. We
  // distinguish between quads and force elements using the element index.

  // Note that we need to provide a force vector (-1, -1) to pull in direction
  // (1, 1).

  for (const auto &element : mesh.localElements()) {
    switch (element.index()) {
    case 0: {
      element_assembler.emplaceElement(
          element, model,
          Integrator(Embedding(Embedding::Collection<Embedding::PhysicalPoint>{{
              vertex_positions.at(connectivity.at(element.index()).at(0)),
              vertex_positions.at(connectivity.at(element.index()).at(1)),
              vertex_positions.at(connectivity.at(element.index()).at(2)),
              vertex_positions.at(connectivity.at(element.index()).at(3)),
          }})));
      break;
    }
    case 1: {
      force_assembler.emplaceElement(element, ForceElement::Force{{-1., -1.}});
      break;
    }
    }
  }

  // We need to create a solver. We do not use the time, so we can set it to
  // zero.

  const auto solver = Solver(&mesh);
  const auto time = Element::Time{0.};

  // Before we can produce meaningful results, we need to specify boundary
  // conditions. Let's fix the nodes at x=0 and y=0.

  std::vector<BoundaryCondition> boundary_conditions;
  for (const auto &vertex : mesh.localVertices()) {
    switch (vertex.index()) {
    case 0:
    case 1:
    case 3: {
      // The displacement in x direction is zero.
      boundary_conditions.push_back({vertex, 0, 0.});
      // The displacement in y direction is zero.
      boundary_conditions.push_back({vertex, 1, 0.});
      break;
    }
    }
  }

  // We are ready to minimize the energy.

  auto result = solver.computeSolution(
      boundary_conditions, Vector::fromGlobalMesh(mesh), time, &assembler);
  cpppetsc::setName("result", &result);

  using DataSource = std::function<void(Mesh::size_type, Mesh::value_type *)>;
  const auto coordinates = cpppetsc::createVectorFromSource(
      mesh, dimension,
      DataSource([&](const Mesh::size_type index, Mesh::value_type *const out) {
        const auto &position = vertex_positions.at(index);
        std::copy(position.begin(), position.end(), out);
      }));

  // Finally we write the results to a file.

  using Viewer = cpppetsc::Viewer<Policy>;
  auto viewer =
      Viewer::fromHdf5FilePath("force_element.ae108", Viewer::Mode::write);
  cpppetsc::writeToViewer(mesh, coordinates, &viewer);
  cpppetsc::writeToViewer(result, &viewer);
}
