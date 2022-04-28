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

add_library(${PROJECT_NAME}-cmdline
            ./CommandLineOptionParser.cc
)

add_library(${PROJECT_NAME}::cmdline ALIAS ${PROJECT_NAME}-cmdline)

target_include_directories(${PROJECT_NAME}-cmdline PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

find_package(Boost 1.67 REQUIRED COMPONENTS program_options)
target_link_libraries(${PROJECT_NAME}-cmdline PUBLIC Boost::program_options)

target_compile_features(${PROJECT_NAME}-cmdline PUBLIC cxx_std_17)

ae108_install_library(cmdline)
