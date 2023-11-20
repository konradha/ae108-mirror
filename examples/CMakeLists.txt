# © 2020 ETH Zurich, Mechanics and Materials Lab
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

project(ae108-examples LANGUAGES C CXX)

if (NOT TARGET ae108::elements)
    find_package(ae108 REQUIRED)
endif()

add_executable(${PROJECT_NAME}-Cmdline Cmdline.cc)
target_link_libraries(${PROJECT_NAME}-Cmdline
                      PRIVATE ae108::cmdline
)

add_executable(${PROJECT_NAME}-MeshGeneration MeshGeneration.cc)
target_link_libraries(${PROJECT_NAME}-MeshGeneration
                      PRIVATE ae108::solve
                      PRIVATE ae108::assembly
                      PRIVATE ae108::elements
)

add_executable(${PROJECT_NAME}-Basic Basic.cc)
target_link_libraries(${PROJECT_NAME}-Basic
                      PRIVATE ae108::solve
                      PRIVATE ae108::assembly
                      PRIVATE ae108::elements
)

add_executable(${PROJECT_NAME}-Output Output.cc)
target_link_libraries(${PROJECT_NAME}-Output
                      PRIVATE ae108::cpppetsc
)

add_executable(${PROJECT_NAME}-Input Input.cc)
target_link_libraries(${PROJECT_NAME}-Input
                      PRIVATE ae108::cpppetsc
)

add_executable(${PROJECT_NAME}-QuadraticTriangles QuadraticTriangles.cc)
target_link_libraries(${PROJECT_NAME}-QuadraticTriangles
                      PRIVATE ae108::solve
                      PRIVATE ae108::assembly
                      PRIVATE ae108::elements
)

add_executable(${PROJECT_NAME}-QuadraticTetrahedra QuadraticTetrahedra.cc)
target_link_libraries(${PROJECT_NAME}-QuadraticTetrahedra
                      PRIVATE ae108::solve
                      PRIVATE ae108::assembly
                      PRIVATE ae108::elements
)

add_executable(${PROJECT_NAME}-SurfaceForce SurfaceForce.cc)
target_link_libraries(${PROJECT_NAME}-SurfaceForce
                      PRIVATE ae108::solve
                      PRIVATE ae108::assembly
                      PRIVATE ae108::elements
)

add_executable(${PROJECT_NAME}-PeriodicBC PeriodicBC.cc)
target_link_libraries(${PROJECT_NAME}-PeriodicBC
                      PRIVATE ae108::solve
                      PRIVATE ae108::assembly
                      PRIVATE ae108::elements
)

add_executable(${PROJECT_NAME}-GeneralizedNonlinearSolver GeneralizedNonlinearSolver.cc)
target_link_libraries(${PROJECT_NAME}-GeneralizedNonlinearSolver
                      PRIVATE ae108::solve
                      PRIVATE ae108::assembly
                      PRIVATE ae108::elements
)

add_executable(${PROJECT_NAME}-TimoshenkoBeam TimoshenkoBeam.cc)
target_link_libraries(${PROJECT_NAME}-TimoshenkoBeam
                      PRIVATE ae108::solve
                      PRIVATE ae108::assembly
                      PRIVATE ae108::elements
)

add_executable(${PROJECT_NAME}-CuboidMesh CuboidMesh.cc)
target_link_libraries(${PROJECT_NAME}-CuboidMesh
                      PRIVATE ae108::cmdline
                      PRIVATE ae108::solve
                      PRIVATE ae108::assembly
                      PRIVATE ae108::elements
)

add_executable(${PROJECT_NAME}-ForceElement ForceElement.cc)
target_link_libraries(${PROJECT_NAME}-ForceElement
                      PRIVATE ae108::cmdline
                      PRIVATE ae108::solve
                      PRIVATE ae108::assembly
                      PRIVATE ae108::elements
)

add_executable(${PROJECT_NAME}-EigenfrequencyAnalysis EigenfrequencyAnalysis.cc)
target_link_libraries(${PROJECT_NAME}-EigenfrequencyAnalysis
                      PRIVATE ae108::cppslepc
                      PRIVATE ae108::assembly
                      PRIVATE ae108::elements
)

if(${AE108_PETSC_COMPLEX})
    add_executable(${PROJECT_NAME}-BlochWaveAnalysis BlochWaveAnalysis.cc)
    target_link_libraries(${PROJECT_NAME}-BlochWaveAnalysis
                          PRIVATE ae108::solve
                          PRIVATE ae108::assembly
                          PRIVATE ae108::elements
                          PRIVATE ae108::cppslepc
    ) 
endif()

add_executable(${PROJECT_NAME}-DynamicSolver DynamicSolver.cc)
target_link_libraries(${PROJECT_NAME}-DynamicSolver
                      PRIVATE ae108::solve
                      PRIVATE ae108::assembly
                      PRIVATE ae108::elements
)

add_executable(${PROJECT_NAME}-CantileverBeam CantileverBeam.cc)
target_link_libraries(${PROJECT_NAME}-CantileverBeam
                      PRIVATE ae108::solve
                      PRIVATE ae108::assembly
                      PRIVATE ae108::elements
)

add_executable(${PROJECT_NAME}-Bar Bar.cc)
target_link_libraries(${PROJECT_NAME}-Bar
                      PRIVATE ae108::solve
                      PRIVATE ae108::assembly
                      PRIVATE ae108::elements
)

if (NOT ${AE108_PETSC_COMPLEX})
    add_executable(${PROJECT_NAME}-NonlinearBC NonlinearBC.cc)
    target_link_libraries(${PROJECT_NAME}-NonlinearBC
                          PRIVATE ae108::solve
                          PRIVATE ae108::assembly
                          PRIVATE ae108::elements
    )
endif()

if (${USE_MESHING})
    add_subdirectory(meshing)
endif()
