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
