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

function(ae108_install_library AE108_LIBRARY)
    include(GNUInstallDirs)
    install(DIRECTORY include/${PROJECT_NAME}
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    )
    set_target_properties(${PROJECT_NAME}-${AE108_LIBRARY}
                          PROPERTIES EXPORT_NAME ${AE108_LIBRARY}
    )
    install(TARGETS ${PROJECT_NAME}-${AE108_LIBRARY}
            EXPORT ${PROJECT_NAME}-${AE108_LIBRARY}-export
            ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
            INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    )
    install(EXPORT ${PROJECT_NAME}-${AE108_LIBRARY}-export
            DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
            NAMESPACE ${PROJECT_NAME}::
            FILE ${PROJECT_NAME}-${AE108_LIBRARY}-export.cmake
    )
endfunction()
