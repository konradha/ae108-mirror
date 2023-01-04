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

add_executable(${PROJECT_NAME}-SolveTests
               ./AffineTransform_Test.cc
               ./DynamicSolver_Test.cc
               ./GeneralizedNonlinearSolver_Test.cc
               ./GeneralizedTAOSolver_Test.cc
               ./InconsistentBoundaryConditionsException_Test.cc
               ./InvalidVertexException_Test.cc
               ./LeastSquaresSolver_Test.cc
               ./NonlinearSolver_Test.cc
               ./TAOSolver_Test.cc
               ./TransformingSolver_Test.cc
               ./boundaryConditionsToEquations_Test.cc
               ./boundaryConditionsToTransform_Test.cc
)

target_include_directories(${PROJECT_NAME}-SolveTests PUBLIC include)

target_link_libraries(${PROJECT_NAME}-SolveTests
                      PRIVATE ${PROJECT_NAME}::solve
                      PRIVATE ${PROJECT_NAME}::cppptest
)

a108_add_cppptest(${PROJECT_NAME}-SolveTests)
