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
