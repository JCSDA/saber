# (C) Copyright 2017-2020 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set(_subdirs bump external gaugrid util)

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
target_link_libraries( ${PROJECT_NAME} PUBLIC eckit )
target_link_libraries( ${PROJECT_NAME} PUBLIC fckit )
if( atlas_FOUND )
    target_link_libraries( ${PROJECT_NAME} PUBLIC atlas_f )
    target_compile_definitions( ${PROJECT_NAME} PRIVATE ATLASIFIED=1)
endif()
if( metis_FOUND )
    target_link_libraries( ${PROJECT_NAME} PUBLIC metis metis_interface )
    target_compile_definitions( ${PROJECT_NAME} PRIVATE SABER_USE_METIS=1 )
endif()

#Fortran module output directory for build and install interfaces
set(MODULE_DIR module/${PROJECT_NAME}/${CMAKE_Fortran_COMPILER_ID}/${CMAKE_Fortran_COMPILER_VERSION})
set_target_properties(${PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR})
install(DIRECTORY ${CMAKE_BINARY_DIR}/${MODULE_DIR}/ DESTINATION ${MODULE_DIR})
target_include_directories(${PROJECT_NAME} INTERFACE
                                $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${MODULE_DIR}>
                                $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/${MODULE_DIR}>)
