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

function(a108_add_cppptest TEST_TARGET)
    add_test(NAME ${TEST_TARGET} COMMAND $<TARGET_FILE:${TEST_TARGET}>)

    find_package(MPI 3.1)
    if (MPI_CXX_FOUND)
        add_test(NAME ${TEST_TARGET}_mpi
                COMMAND ${MPIEXEC_EXECUTABLE}
                        ${MPIEXEC_NUMPROC_FLAG} 2
                        ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${TEST_TARGET}> ${MPIEXEC_POSTFLAGS}
        )
    endif()
endfunction()
