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

add_library(${PROJECT_NAME}-elements
            AutomaticForcesTrait.cc
            AutomaticStiffnessMatrixTrait.cc
            Bar.cc
            ComputeConsistentMassMatrixTrait.cc
            ComputeEnergyTrait.cc
            ComputeForcesTrait.cc
            ComputeLumpedMassMatrixTrait.cc
            ComputeStiffnessMatrixTrait.cc
            ComputeMassMatrixTrait.cc
            CoreElement.cc
            TwoNodeCorotationalBeamElement.cc
            ElementBase.cc
            ElementWithMass.cc
            ForceElement.cc
            Minimal.cc
            TimoshenkoBeamElement.cc
            automatic_forces.cc
            automatic_stiffness_matrix.cc
            compute_consistent_mass_matrix.cc
            compute_energy.cc
            compute_forces.cc
            compute_lumped_mass_matrix.cc
            compute_stiffness_matrix.cc
            embedding/AutomaticJacobianTrait.cc
            embedding/EmbedTrait.cc
            embedding/EmbeddingBase.cc
            embedding/IsoparametricEmbedding.cc
            embedding/JacobianTrait.cc
            embedding/automatic_jacobian.cc
            embedding/compute_jacobian.cc
            embedding/embed_point.cc
            integrator/ComputeVolumeTrait.cc
            integrator/IntegrateShapeTrait.cc
            integrator/IntegrateTrait.cc
            integrator/IntegratorBase.cc
            integrator/IsoparametricIntegrator.cc
            integrator/compute_volume.cc
            integrator/integrate.cc
            integrator/integrate_shape.cc
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
            mesh/generate_cuboid_mesh.cc
            mesh/generate_quadratic_tetrahedron_mesh.cc
            mesh/generate_quadratic_triangle_mesh.cc
            mesh/generate_tetrahedron_mesh.cc
            mesh/generate_triangle_mesh.cc
            mesh/refine_segment_mesh.cc
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
            shape/ShapeBase.cc
            shape/Tet10.cc
            shape/Tet4.cc
            shape/Tri3.cc
            shape/Tri6.cc
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
                        PUBLIC cxx_std_17
)

find_package(Eigen3 3.4 CONFIG REQUIRED)
find_package(Boost 1.67 REQUIRED)

target_link_libraries(${PROJECT_NAME}-elements
                      PUBLIC Eigen3::Eigen
                      PUBLIC Boost::boost
)

add_library(${PROJECT_NAME}::elements ALIAS ${PROJECT_NAME}-elements)

ae108_install_library(elements)
