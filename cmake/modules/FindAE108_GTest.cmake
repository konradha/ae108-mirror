# © 2021 ETH Zurich, Mechanics and Materials Lab
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

find_package(GTest 1.8.1 CONFIG)
if (NOT GTest_FOUND)
    set(FETCHCONTENT_SOURCE_DIR_AE108_GOOGLETEST "/usr/src/googletest")
    if (IS_DIRECTORY "${FETCHCONTENT_SOURCE_DIR_AE108_GOOGLETEST}")
        message(STATUS
            "GTest was not found, using source at ${FETCHCONTENT_SOURCE_DIR_AE108_GOOGLETEST} instead."
        )

        include(FetchContent)
        FetchContent_Declare("ae108_googletest")
        FetchContent_GetProperties(
            ae108_googletest
            SOURCE_DIR ae108_googletest_SOURCE_DIR
            BINARY_DIR ae108_googletest_BINARY_DIR
            POPULATED ae108_googletest_POPULATED
        )
        if (NOT ae108_googletest_POPULATED)
            FetchContent_Populate("ae108_googletest")
            add_subdirectory(
                "${ae108_googletest_SOURCE_DIR}"
                "${ae108_googletest_BINARY_DIR}"
                EXCLUDE_FROM_ALL
            )
        endif()

        foreach(LIBRARY_NAME gtest gmock gtest_main gmock_main)
            add_library("GTest::${LIBRARY_NAME}" ALIAS "${LIBRARY_NAME}")
        endforeach()

        set(GTest_FOUND TRUE)
    endif()
endif()

find_package_handle_standard_args(AE108_GTest
    REQUIRED_VARS GTest_FOUND
)