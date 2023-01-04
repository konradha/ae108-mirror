# © 2020 ETH Zurich, Mechanics and Materials Lab
# © 2020 California Institute of Technology
#
# This file is part of ae108.
#
# ae108 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or any
# later version.
#
# ae108 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ae108. If not, see <https://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

add_library(${PROJECT_NAME}-cpppetsc
            ./Context.cc
            ./Exception.cc
            ./GeneralizedMeshBoundaryConditition.cc
            ./InvalidParametersException.cc
            ./LeastSquaresSolver.cc
            ./LinearSolver.cc
            ./LinearSolverDivergedException.cc
            ./LocalVertexView.cc
            ./Matrix.cc
            ./Matrix_fwd.cc
            ./Mesh.cc
            ./MeshBoundaryCondition.cc
            ./MeshDataProvider.cc
            ./Mesh_fwd.cc
            ./NonlinearSolver.cc
            ./NonlinearSolverDivergedException.cc
            ./ParallelComputePolicy.cc
            ./ParallelComputePolicy_fwd.cc
            ./SequentialComputePolicy.cc
            ./SequentialComputePolicy_fwd.cc
            ./SharedEntity.cc
            ./TAOSolver.cc
            ./TAOSolverDivergedException.cc
            ./TaggedEntity.cc
            ./TaggedVector.cc
            ./UniqueEntity.cc
            ./Vector.cc
            ./Vector_fwd.cc
            ./Viewer.cc
            ./asHermitianTransposedMatrix.cc
            ./asInverseMatrix.cc
            ./asSchurComplement.cc
            ./asThAT.cc
            ./asTransposedMatrix.cc
            ./clone.cc
            ./computeElementsOfMatrix.cc
            ./computeSchurComplement.cc
            ./copy.cc
            ./createLhsTransform.cc
            ./createRhsTransform.cc
            ./createTransformInput.cc
            ./createTransformOutput.cc
            ./createVectorFromSource.cc
            ./getName.cc
            ./multiply.cc
            ./nestMatrices.cc
            ./nestVectors.cc
            ./readFromViewer.cc
            ./scalarProduct.cc
            ./setName.cc
            ./vertexDataOffsets.cc
            ./writeToViewer.cc
)

add_library(${PROJECT_NAME}::cpppetsc ALIAS ${PROJECT_NAME}-cpppetsc)

target_include_directories(${PROJECT_NAME}-cpppetsc PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

find_package(AE108_PETSc 3.15 MODULE REQUIRED)
find_package(range-v3 0.11.0 CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME}-cpppetsc
                      PUBLIC ${PROJECT_NAME}::external::petsc
                      PUBLIC range-v3::range-v3
)

target_compile_features(${PROJECT_NAME}-cpppetsc PUBLIC cxx_std_17)

ae108_install_library(cpppetsc)
