# (C) Copyright 2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list( APPEND saber_quench_src_files
Fields.h
Fields.cc
Geometry.h
Geometry.cc
Increment.h
Increment.cc
State.h
State.cc
)

ecbuild_add_library( TARGET   quench
                     SOURCES  ${saber_quench_src_files}
                     PUBLIC_LIBS  NetCDF::NetCDF_CXX oops
                     HEADER_DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/quench
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE CXX )

#Configure include directory layout for build-tree to match install-tree
set(QUENCH_BUILD_DIR_INCLUDE_PATH ${CMAKE_BINARY_DIR}/${PROJECT_NAME}/include)
add_custom_target(quench_headers ALL COMMAND ${CMAKE_COMMAND} -E make_directory "${QUENCH_BUILD_DIR_INCLUDE_PATH}/saber"
                                     COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR} ${QUENCH_BUILD_DIR_INCLUDE_PATH}/saber/quench)
add_dependencies(quench quench_headers)
target_include_directories(quench PUBLIC $<BUILD_INTERFACE:${QUENCH_BUILD_DIR_INCLUDE_PATH}>
                                         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
                                         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../..>)

ecbuild_add_test( TARGET quench_coding_norms
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/oops_cpplint.py
                  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR} )