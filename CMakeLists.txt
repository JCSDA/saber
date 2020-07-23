# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# SABER
################################################################################

cmake_minimum_required( VERSION 3.12 )

project( saber LANGUAGES CXX Fortran )

## Ecbuild integration
find_package( ecbuild QUIET )
include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_declare_project()
list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
set( CMAKE_DIRECTORY_LABELS ${PROJECT_NAME} )

## Configuration options
set( ENABLE_SABER_DOC "Build SABER/BUMP Documentation" OFF )
include( ${PROJECT_NAME}_compiler_flags )

## Dependencies
# Required
find_package( OpenMP REQUIRED COMPONENTS Fortran )
find_package( MPI REQUIRED COMPONENTS Fortran )
find_package( NetCDF REQUIRED COMPONENTS Fortran )
find_package( eckit REQUIRED COMPONENTS MPI )
find_package( fckit REQUIRED )
find_package( atlas REQUIRED COMPONENTS OMP_Fortran )
# Optionals
find_package( oops QUIET OPTIONAL_COMPONENTS qg )

## Sources
add_subdirectory( src/saber )
add_subdirectory( test )
add_subdirectory( tools )
add_subdirectory( doc )

## Global tests
ecbuild_add_test( TARGET saber_coding_norms
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/${PROJECT_NAME}_cpplint.py
                  ARGS --quiet --recursive ${CMAKE_SOURCE_DIR} )

## Package Config
ecbuild_install_project( NAME ${PROJECT_NAME} )

ecbuild_print_summary()
