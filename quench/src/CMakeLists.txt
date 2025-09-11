# (C) Copyright 2022 UCAR.
# (C) Copyright 2025 Meteorologisk Institutt
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

message( STATUS "Found quench I/O formats:" )
message( STATUS "- Default (NetCDF)" )
message( STATUS "- GMSH" )

# Common sources, not obs-related
list( APPEND quench_src_files

Covariance.h
Fields.cc
Fields.h
Geometry.cc
Geometry.h
Increment.cc
Increment.h
Interpolation.cc
Interpolation.h
LinearVariableChange.h
LinearVariableChange.cc
LinearVariableChangeParameters.h
ModelData.h
State.cc
State.h
VariableChange.h
VariableChange.cc

FieldsIO/FieldsIOBase.cc
FieldsIO/FieldsIOBase.h
FieldsIO/FieldsIODefault.cc
FieldsIO/FieldsIODefault.h
FieldsIO/FieldsIOGmsh.cc
FieldsIO/FieldsIOGmsh.h

)

ecbuild_add_library( TARGET   quench
                     SOURCES  ${quench_src_files}
                     PUBLIC_LIBS  oops saber
                     HEADER_DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/quench
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE CXX )

target_link_libraries( quench PUBLIC NetCDF::NetCDF_Fortran NetCDF::NetCDF_C )
if( eccodes_FOUND )
    target_link_libraries( quench PUBLIC eccodes )
endif()

#Configure include directory layout for build-tree to match install-tree
set(QUENCH_BUILD_DIR_INCLUDE_PATH ${CMAKE_BINARY_DIR}/${PROJECT_NAME}/include)
add_custom_target(quench_headers ALL COMMAND ${CMAKE_COMMAND} -E make_directory "${QUENCH_BUILD_DIR_INCLUDE_PATH}/saber"
                                     COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR} ${QUENCH_BUILD_DIR_INCLUDE_PATH}/saber/quench)
add_dependencies(quench quench_headers)
target_include_directories(quench PUBLIC $<BUILD_INTERFACE:${QUENCH_BUILD_DIR_INCLUDE_PATH}>
                                         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
                                         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../..>)
