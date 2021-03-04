# © 2020 ETH Zurich, Mechanics and Materials Lab
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

add_library(${PROJECT_NAME}-elements
            AutomaticForcesTrait.cc
            AutomaticStiffnessMatrixTrait.cc
            ComputeEnergyTrait.cc
            ComputeForcesTrait.cc
            ComputeStiffnessMatrixTrait.cc
            CoreElement.cc
            ElementBase.cc
            ForceElement.cc
            Minimal.cc
            automatic_forces.cc
            automatic_stiffness_matrix.cc
            compute_energy.cc
            compute_forces.cc
            compute_stiffness_matrix.cc
            embedding/AutomaticJacobianTrait.cc
            embedding/EmbedTrait.cc
            embedding/EmbeddingBase.cc
            embedding/IsoparametricEmbedding.cc
            embedding/JacobianTrait.cc
            embedding/automatic_jacobian.cc
            embedding/compute_jacobian.cc
            embedding/embed_point.cc
            integrator/IntegrateTrait.cc
            integrator/IntegratorBase.cc
            integrator/IsoparametricIntegrator.cc
            integrator/integrate.cc
            materialmodels/AutomaticStressTrait.cc
            materialmodels/AutomaticTangentMatrixTrait.cc
            materialmodels/ComputeEnergyTrait.cc
            materialmodels/ComputeStrainTrait.cc
            materialmodels/ComputeStressTrait.cc
            materialmodels/ComputeTangentMatrixTrait.cc
            materialmodels/Hookean.cc
            materialmodels/MaterialModelBase.cc
            materialmodels/Minimal.cc
            materialmodels/automatic_stress.cc
            materialmodels/automatic_tangent_matrix.cc
            materialmodels/compute_energy.cc
            materialmodels/compute_strain.cc
            materialmodels/compute_stress.cc
            materialmodels/compute_tangent_matrix.cc
            mesh/Mesh.cc
            mesh/generate_quadratic_triangle_mesh.cc
            mesh/generate_tetrahedron_mesh.cc
            mesh/generate_triangle_mesh.cc
            quadrature/IntegrateTrait.cc
            quadrature/Quadrature.cc
            quadrature/QuadratureBase.cc
            quadrature/integrate.cc
            shape/AutomaticGradientTrait.cc
            shape/GradientTrait.cc
            shape/Hexa8.cc
            shape/PointTrait.cc
            shape/Quad4.cc
            shape/Seg2.cc
            shape/Tet10.cc
            shape/Tet4.cc
            shape/Tri3.cc
            shape/Tri6.cc
            shape/ShapeBase.cc
            shape/ValueTrait.cc
            shape/automatic_gradients.cc
            shape/compute_gradients.cc
            shape/compute_values.cc
            shape/get_points.cc
            tensor/Tensor.cc
            tensor/as_matrix_of_columns.cc
            tensor/as_matrix_of_rows.cc
            tensor/as_two_tensor.cc
            tensor/as_vector.cc
            tensor/differentiate.cc
)

target_include_directories(${PROJECT_NAME}-elements PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_compile_features(${PROJECT_NAME}-elements
                        PUBLIC cxx_std_11
)

find_package(Eigen3 3.3 CONFIG REQUIRED)
find_package(Boost 1.67 REQUIRED)

target_link_libraries(${PROJECT_NAME}-elements
                      PUBLIC Eigen3::Eigen
                      PUBLIC Boost::boost
)

add_library(${PROJECT_NAME}::elements ALIAS ${PROJECT_NAME}-elements)

ae108_install_library(elements)
