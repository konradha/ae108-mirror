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

add_library(${PROJECT_NAME}-cppptest STATIC
            ./Matchers.cc
            ./isLocal.cc
            ./main.cc
)

add_library(${PROJECT_NAME}::cppptest ALIAS ${PROJECT_NAME}-cppptest)

target_include_directories(${PROJECT_NAME}-cppptest PUBLIC include)

find_package(GTest 1.8.1 CONFIG REQUIRED)
find_package(AE108_PETSc 3.12 MODULE REQUIRED)
target_link_libraries(${PROJECT_NAME}-cppptest
                      PRIVATE ${PROJECT_NAME}::external::petsc
                      PUBLIC GTest::gmock
)

target_compile_features(${PROJECT_NAME}-cppptest PUBLIC cxx_std_17)
