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
#include "ae108/elements/CoreElement.h"
#include "ae108/elements/embedding/IsoparametricEmbedding.h"
#include "ae108/elements/integrator/IsoparametricIntegrator.h"
#include "ae108/elements/materialmodels/Hookean.h"
#include "ae108/elements/mesh/generate_quadratic_triangle_mesh.h"
#include "ae108/elements/quadrature/Quadrature.h"
#include "ae108/elements/shape/Tri6.h"
#include "ae108/solve/NonlinearSolver.h"
#include <array>

using namespace ae108;

using Policy = cpppetsc::ParallelComputePolicy;
using Context = cpppetsc::Context<Policy>;
using Mesh = cpppetsc::Mesh<Policy>;
using Vector = cpppetsc::Vector<Policy>;
using Viewer = cpppetsc::Viewer<Policy>;
using BoundaryCondition = cpppetsc::MeshBoundaryCondition<Mesh>;

namespace embedding = ae108::elements::embedding;
namespace integrator = ae108::elements::integrator;
namespace materialmodels = ae108::elements::materialmodels;
namespace quadrature = ae108::elements::quadrature;
namespace shape = ae108::elements::shape;
namespace mesh = ae108::elements::mesh;

// In this example we will simulate a mesh of quadratic triangles that
// that we will generate later.

// Let's specify the parameters of this mesh.

constexpr auto number_of_vertices_per_element = Mesh::size_type{6};
constexpr auto dimension = Mesh::size_type{2};
constexpr auto dof_per_vertex = Mesh::size_type{2};

// Let's specify how the elements behave when deformed.

// First, we select a Hookean (linear elastic) material model in 2D.

using MaterialModel =
    materialmodels::Hookean<dimension, Vector::value_type, Vector::real_type>;

// We'll choose the Tri6 shape functions (6 shape functions) that are defined
// in 2D reference space.

using Shape = shape::Tri6;

// Since also our physical space is 2D we may use the isoparametric embedding.

using Embedding = embedding::IsoparametricEmbedding<Shape>;

// Let's choose a simple quadrature rule in reference space.

using Quadrature = quadrature::Quadrature<quadrature::QuadratureType::Simplex,
                                          dimension, 2 /* order */>;

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

// Our goal is to minimize the energy. This is done by the solver.

using Solver = solve::NonlinearSolver<Assembler>;

int main(int argc, char **argv) {
  // MPI/PETSc/cpppetsc must be initialized before using it.

  const auto context = Context(&argc, &argv);

  // Now we generate a geometry of 4 elements connected in x-direction.
  //  *-----*-----*
  //  |  \  |  \  |
  //  *-----*-----*
  const auto geometry =
      mesh::generate_quadratic_triangle_mesh({{2., 1.}}, {{2, 1}});

  const auto mesh =
      Mesh::fromConnectivity(dimension, geometry.connectivity(),
                             geometry.number_of_positions(), dof_per_vertex);
  auto assembler = Assembler();

  const auto model = MaterialModel(1.0, 0.);

  // Depending on whether we use MPI, our mesh may be distributed and not all
  // elements are present on this computational node.

  // Let's add those elements that are "local" to the assembler.

  for (const auto &element : mesh.localElements()) {
    const auto vertexIndices = element.vertexIndices();
    assembler.emplaceElement(
        element, model,
        Integrator(Embedding(Embedding::Collection<Embedding::PhysicalPoint>{{
            geometry.position_of_vertex(vertexIndices.at(0)),
            geometry.position_of_vertex(vertexIndices.at(1)),
            geometry.position_of_vertex(vertexIndices.at(2)),
            geometry.position_of_vertex(vertexIndices.at(3)),
            geometry.position_of_vertex(vertexIndices.at(4)),
            geometry.position_of_vertex(vertexIndices.at(5)),
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
    const auto position = geometry.position_of_vertex(vertex.index());
    const auto tolerance = .1;
    if (std::abs(position[0] - 0.) < tolerance) {
      // The displacement in x direction is zero.
      boundary_conditions.push_back({vertex, 0, 0.});
      // The displacement in y direction is zero.
      boundary_conditions.push_back({vertex, 1, 0.});
    } else if (std::abs(position[0] - 2.) < tolerance) {
      // The displacement in x direction is .5.
      boundary_conditions.push_back({vertex, 0, .5});
      // The displacement in y direction is 0.
      boundary_conditions.push_back({vertex, 1, 0.});
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
      mesh, dimension,
      DataSource([&](const Mesh::size_type index, Mesh::value_type *const out) {
        const auto &position = geometry.position_of_vertex(index);
        std::copy(position.begin(), position.end(), out);
      }));
  cpppetsc::setName("coordinates", &coordinates);

  // Now we write the mesh to a file.
  auto viewer = Viewer::fromHdf5FilePath("quadratic_triangles.ae108",
                                         Viewer::Mode::write);
  cpppetsc::writeToViewer(mesh, coordinates, &viewer);

  // Let's write the result to the file.
  cpppetsc::setName("result", &result);
  cpppetsc::writeToViewer(result, &viewer);

  fprintf(stderr, "The data has been written to the file "
                  "\"quadratic_triangles.ae108\".\n");
}
