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
#include "ae108/cmdline/CommandLineOptionParser.h"
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
#include "ae108/elements/mesh/generate_cuboid_mesh.h"
#include "ae108/elements/quadrature/Quadrature.h"
#include "ae108/elements/shape/Hexa8.h"
#include "ae108/solve/NonlinearSolver.h"
#include <iostream>

using namespace ae108;

using Policy = cpppetsc::ParallelComputePolicy;
using Context = cpppetsc::Context<Policy>;
using Mesh = cpppetsc::Mesh<Policy>;
using Vector = cpppetsc::Vector<Policy>;
using BoundaryCondition = cpppetsc::MeshBoundaryCondition<Mesh>;
using Viewer = cpppetsc::Viewer<Policy>;

namespace embedding = ae108::elements::embedding;
namespace integrator = ae108::elements::integrator;
namespace materialmodels = ae108::elements::materialmodels;
namespace quadrature = ae108::elements::quadrature;
namespace shape = ae108::elements::shape;
namespace mesh = ae108::elements::mesh;

// In this example we will simulate a mesh of cuboids that we will generate
// later.

// Let's specify the parameters of this mesh.

constexpr auto dimension = Mesh::size_type{3};
constexpr auto dof_per_vertex = Mesh::size_type{3};

// Let's specify how the elements behave when deformed.

// First, we select a Hookean (linear elastic) material model in 3D.

using MaterialModel =
    materialmodels::Hookean<dimension, Vector::value_type, Vector::real_type>;

// We'll choose the Hexa8 shape functions (8 shape functions) that are defined
// in 3D reference space.

using Shape = shape::Hexa8;

// Since also our physical space is 3D we may use the isoparametric embedding.

using Embedding = embedding::IsoparametricEmbedding<Shape>;

// Let's choose a simple quadrature rule in reference space.

using Quadrature = quadrature::Quadrature<quadrature::QuadratureType::Cube,
                                          dimension, 3 /* order */>;

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
  // First we parse the command line flags "granularity", "stdout-output", and
  // "file-output".

  auto granularity = std::size_t{2};
  auto stdout_output = true;
  auto file_output = false;
  cmdline::CommandLineOptionParser{std::cerr}
      .withOption("granularity",
                  "the number of cuboids per coordinates axis (default: 2)",
                  &granularity)
      .withOption("stdout-output", "enable output to stdout (default: true)",
                  &stdout_output)
      .withOption("file-output", "enable output to HDF5 file (default: false)",
                  &file_output)
      .parse(argc, argv);

  // MPI/PETSc/cpppetsc must be initialized before using it.

  const auto context = Context(&argc, &argv);

  // Now we generate a geometry of granularity^3 cuboids.

  const auto geometry = mesh::generate_cuboid_mesh(
      {{1., 1., 1.}}, {{granularity, granularity, granularity}});

  const auto mesh =
      Mesh::fromConnectivity(dimension, geometry.connectivity(),
                             geometry.number_of_positions(), dof_per_vertex);
  auto assembler = Assembler();

  const auto model = MaterialModel(1.0, .1);

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
            geometry.position_of_vertex(vertexIndices.at(6)),
            geometry.position_of_vertex(vertexIndices.at(7)),
        }})));
  }

  // We need to create a solver. We do not use the time, so we can set it to
  // zero.

  const auto solver = Solver(&mesh);
  const auto time = Element::Time{0.};

  // Before we can produce meaningful results, we need to specify boundary
  // conditions. Let's fix the nodes at x=0 and pull on the nodes at x=1.

  std::vector<BoundaryCondition> boundary_conditions;
  for (const auto &vertex : mesh.localVertices()) {
    const auto position = geometry.position_of_vertex(vertex.index());
    const auto tolerance = 1. / static_cast<double>(granularity) / 2.;
    if (std::abs(position[0] - 0.) < tolerance) {
      // The displacement in x direction is zero.
      boundary_conditions.push_back({vertex, 0, 0.});
      // The displacement in y direction is zero.
      boundary_conditions.push_back({vertex, 1, 0.});
      // The displacement in z direction is zero.
      boundary_conditions.push_back({vertex, 2, 0.});
    } else if (std::abs(position[0] - 1.) < tolerance) {
      // The displacement in x direction is .5.
      boundary_conditions.push_back({vertex, 0, .5});
      // The displacement in y direction is 0.
      boundary_conditions.push_back({vertex, 1, 0.});
      // The displacement in z direction is 0.
      boundary_conditions.push_back({vertex, 2, 0.});
    }
  }

  // We are ready to minimize the energy.

  auto result = solver.computeSolution(
      boundary_conditions, Vector::fromGlobalMesh(mesh), time, &assembler);
  cpppetsc::setName("result", &result);

  // As a final step, we print to result to stdout and/or write it to a file.

  if (stdout_output) {
    const auto global_result =
        Mesh::vector_type::fromDistributedInCanonicalOrder(result, mesh);
    if (Policy::isPrimaryRank()) {
      global_result.unwrap().print();
    }
  }

  if (file_output) {
    using DataSource = std::function<void(Mesh::size_type, Mesh::value_type *)>;
    const auto coordinates = cpppetsc::createVectorFromSource(
        mesh, dimension,
        DataSource(
            [&](const Mesh::size_type index, Mesh::value_type *const out) {
              const auto &position = geometry.position_of_vertex(index);
              std::copy(position.begin(), position.end(), out);
            }));

    auto viewer =
        Viewer::fromHdf5FilePath("cuboid_mesh.ae108", Viewer::Mode::write);
    cpppetsc::writeToViewer(mesh, coordinates, &viewer);
    cpppetsc::writeToViewer(result, &viewer);
  }
}
