# (C) Copyright 2017-2020 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Build list of subdirs with files to add
set(_subdirs bump external gaugrid util)
if( oops_FOUND )
    list(APPEND _subdirs oops interpolation)
endif()
foreach( _subdir IN LISTS _subdirs )
    add_subdirectory( ${_subdir} )
    list( TRANSFORM ${_subdir}_src_files PREPEND ${_subdir}/ )
    list( APPEND ${PROJECT_NAME}_src_files ${${_subdir}_src_files} )
endforeach()

include(GNUInstallDirs)
ecbuild_add_library( TARGET ${PROJECT_NAME}
                     SOURCES ${${PROJECT_NAME}_src_files}
                     INSTALL_HEADERS LISTED
                     HEADER_DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
                     LINKER_LANGUAGE CXX )

target_link_libraries( ${PROJECT_NAME} PUBLIC OpenMP::OpenMP_Fortran )
target_link_libraries( ${PROJECT_NAME} PUBLIC NetCDF::NetCDF_Fortran )
target_link_libraries( ${PROJECT_NAME} PUBLIC MPI::MPI_Fortran )
target_link_libraries( ${PROJECT_NAME} PUBLIC ${LAPACK_LIBRARIES} )
target_link_libraries( ${PROJECT_NAME} PUBLIC eckit )
target_link_libraries( ${PROJECT_NAME} PUBLIC fckit )
target_link_libraries( ${PROJECT_NAME} PUBLIC atlas_f )
if( oops_FOUND )
    target_link_libraries( ${PROJECT_NAME} PUBLIC ${oops_LIBRARIES} )
endif()

# Enable including of local header Fortran file templates
target_include_directories( ${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..> )

# Process fypp files
fckit_target_preprocess_fypp( ${PROJECT_NAME}
                              DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/generics.fypp )

# Fortran module output directory for build and install interfaces
set( MODULE_DIR module/${PROJECT_NAME}/${CMAKE_Fortran_COMPILER_ID}/${CMAKE_Fortran_COMPILER_VERSION} )
set_target_properties( ${PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR} )
install( DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR}/ DESTINATION ${MODULE_DIR} )
target_include_directories( ${PROJECT_NAME} INTERFACE
                            $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${MODULE_DIR}>
                            $<INSTALL_INTERFACE:${MODULE_DIR}> )
