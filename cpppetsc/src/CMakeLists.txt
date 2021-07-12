# © 2020 ETH Zurich, Mechanics and Materials Lab
# © 2020 California Institute of Technology
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

add_library(${PROJECT_NAME}-cpppetsc
            ./Context.cc
            ./Exception.cc
            ./GeneralizedMeshBoundaryConditition.cc
            ./IndexBasedAccessIterator.cc
            ./InvalidParametersException.cc
            ./IteratorRange.cc
            ./LinearSolver.cc
            ./LinearSolverDivergedException.cc
            ./LocalVertexView.cc
            ./Matrix.cc
            ./Matrix_fwd.cc
            ./Mesh.cc
            ./MeshBoundaryCondition.cc
            ./Mesh_fwd.cc
            ./NonlinearSolver.cc
            ./NonlinearSolverDivergedException.cc
            ./ParallelComputePolicy.cc
            ./ParallelComputePolicy_fwd.cc
            ./SequentialComputePolicy.cc
            ./SequentialComputePolicy_fwd.cc
            ./TAOSolver.cc
            ./TAOSolverDivergedException.cc
            ./TaggedEntity.cc
            ./TaggedVector.cc
            ./UniqueEntity.cc
            ./Vector.cc
            ./Vector_fwd.cc
            ./Viewer.cc
            ./asInverseMatrix.cc
            ./asSchurComplement.cc
            ./asTransformedMatrix.cc
            ./clone.cc
            ./copy.cc
            ./createRhsTransform.cc
            ./createLhsTransform.cc
            ./createTransformInput.cc
            ./createTransformOutput.cc
            ./createVectorFromSource.cc
            ./getName.cc
            ./multiply.cc
            ./readFromViewer.cc
            ./setName.cc
            ./vertexDataOffsets.cc
            ./writeToViewer.cc
)

add_library(${PROJECT_NAME}::cpppetsc ALIAS ${PROJECT_NAME}-cpppetsc)

target_include_directories(${PROJECT_NAME}-cpppetsc PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

find_package(AE108_PETSc 3.12 MODULE REQUIRED)
find_package(range-v3 0.10.0 CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME}-cpppetsc
                      PUBLIC ${PROJECT_NAME}::external::petsc
                      PUBLIC range-v3::range-v3
)

target_compile_features(${PROJECT_NAME}-cpppetsc PUBLIC cxx_std_11)

ae108_install_library(cpppetsc)
