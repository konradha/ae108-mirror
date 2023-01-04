# © 2021 ETH Zurich, Mechanics and Materials Lab
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

add_executable(${PROJECT_NAME}-CppSlepcTests
        EigenvalueProblemSolverDivergedException_Test.cc
        InvalidEigenvalueIndexException_Test.cc
        InvalidProblemTypeException_Test.cc
        LinearEigenvalueProblemSolver_Test.cc
        NoOperatorsSetException_Test.cc
        computeEigenvalues_Test.cc
        computeSmallestEigenvalues_Test.cc
        main.cc
)

add_test(NAME ${PROJECT_NAME}-CppSlepcTests
         COMMAND $<TARGET_FILE:${PROJECT_NAME}-CppSlepcTests>)

find_package(GTest 1.8.1 CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME}-CppSlepcTests
                      PRIVATE ${PROJECT_NAME}::cppslepc
                      PRIVATE GTest::gmock
                      PRIVATE ${PROJECT_NAME}::cppptest
)
