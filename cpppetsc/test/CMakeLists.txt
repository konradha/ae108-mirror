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

add_executable(${PROJECT_NAME}-CppPetscTests
               ./InvalidParametersException_Test.cc
               ./LeastSquaresSolver_Test.cc
               ./LinearSolverDivergedException_Test.cc
               ./LinearSolver_Test.cc
               ./Matrix_Test.cc
               ./Mesh_Test.cc
               ./NonlinearSolverDivergedException_Test.cc
               ./NonlinearSolver_Test.cc
               ./SharedEntity_Test.cc
               ./TAOSolverDivergedException_Test.cc
               ./TAOSolver_Test.cc
               ./TaggedEntity_Test.cc
               ./Vector_Test.cc
               ./asHermitianTransposedMatrix_Test.cc
               ./asInverseMatrix_Test.cc
               ./asSchurComplement_Test.cc
               ./asThAT_Test.cc
               ./asTransposedMatrix_Test.cc
               ./clone_Test.cc
               ./computeElementsOfMatrix_Test.cc
               ./computeSchurComplement_Test.cc
               ./copy_Test.cc
               ./createLhsTransform_Test.cc
               ./createRhsTransform_Test.cc
               ./createTransformInput_Test.cc
               ./createTransformOutput_Test.cc
               ./createVectorFromSource_Test.cc
               ./multiply_Test.cc
               ./nestMatrices_Test.cc
               ./nestVectors_Test.cc
               ./scalarProduct_Test.cc
               ./setName_Test.cc
               ./vertexDataOffsets_Test.cc
)

target_include_directories(${PROJECT_NAME}-CppPetscTests PRIVATE include)

target_link_libraries(${PROJECT_NAME}-CppPetscTests
                      PRIVATE ${PROJECT_NAME}::cpppetsc
                      PRIVATE ${PROJECT_NAME}::cppptest
)

a108_add_cppptest(${PROJECT_NAME}-CppPetscTests)
