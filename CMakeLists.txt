# (C) Copyright 2017-2020 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# SABER
################################################################################

cmake_minimum_required( VERSION 3.12 )

project( saber VERSION 1.3.0 LANGUAGES CXX Fortran )

option(OPENMP "Build saber with OpenMP support" ON)

## Ecbuild integration
find_package( ecbuild QUIET )
include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_declare_project()
list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
set( CMAKE_DIRECTORY_LABELS ${PROJECT_NAME} )

## Configuration options
option( ENABLE_SABER_DOC "Build SABER documentation" OFF )
option( ENABLE_MKL "Use MKL for LAPACK implementation (if available)" ON )
include( ${PROJECT_NAME}_compiler_flags )

## Dependencies
# Required
find_package( jedicmake REQUIRED )  # Prefer find modules from jedi-cmake
if(OPENMP)
  find_package( OpenMP REQUIRED COMPONENTS Fortran )
endif()
find_package( MPI REQUIRED COMPONENTS Fortran )
find_package( NetCDF REQUIRED COMPONENTS C Fortran )
find_package( eckit 1.11.6 REQUIRED COMPONENTS MPI )
find_package( fckit 0.7.0  REQUIRED )
find_package( atlas 0.27.0 REQUIRED COMPONENTS OMP_Fortran )

if( ENABLE_MKL )
    find_package( MKL )
endif()
if( MKL_FOUND )
    set( LAPACK_LIBRARIES ${MKL_LIBRARIES} )
else()
    find_package( LAPACK REQUIRED )
endif()

# Optionals
find_package( oops QUIET OPTIONAL_COMPONENTS qg )
find_package( gsibec 1.0.4 QUIET )
find_package( vader QUIET )

# Optional SABER blocks
if( gsibec_FOUND )
    find_package( sp QUIET )
    message( STATUS "SABER block GSI is enabled" )
else()
    message( STATUS "SABER block GSI is not enabled" )
endif()
# TODO(Benjamin): remove this line once ATLAS is upgraded to 0.30.0 everywhere
set( ENABLE_SPECTRALB 0 )
if( atlas_VERSION GREATER_EQUAL 0.30.0 )
    message( STATUS "ATLAS version higher or equal to 0.30.0" )
    if( atlas_TRANS_FOUND )
        message( STATUS "SABER block SPECTRALB is enabled" )
        add_definitions(-Datlas_TRANS_FOUND=1)
        set( ENABLE_SPECTRALB 1 )
    else()
        message( STATUS "SABER block SPECTRALB is not enabled" )
        add_definitions(-Datlas_TRANS_FOUND=0)
    endif()
else()
    message( STATUS "ATLAS version lower than 0.30.0" )
    message( STATUS "SABER block SPECTRALB is not enabled" )
    add_definitions(-Datlas_TRANS_FOUND=0)
endif()
if( vader_FOUND )
    if ( ENABLE_VADER_MO )
      message( STATUS "SABER blocks VADER are enabled (with MO-specific blocks)" )
    else()
      message( STATUS "SABER blocks VADER are enabled (without MO-specific blocks)" )
    endif()
else()
    message( STATUS "SABER blocks VADER are not enabled" )
endif()

## SABER instrumentation
set( ENABLE_SABER_INSTRUMENTATION 0 )
if( DEFINED ENV{ENABLE_SABER_INSTRUMENTATION} )
    set( ENABLE_SABER_INSTRUMENTATION $ENV{ENABLE_SABER_INSTRUMENTATION} )
endif()
if( ENABLE_SABER_INSTRUMENTATION )
    message( STATUS "SABER instrumentation is enabled" )
    add_definitions(-DENABLE_SABER_INSTRUMENTATION=1)
else()
    message( STATUS "SABER instrumentation is not enabled" )
    add_definitions(-DENABLE_SABER_INSTRUMENTATION=0)
endif()

## Sources
add_subdirectory( src/saber )
add_subdirectory( tools ) # tools is needed for test
add_subdirectory( test )
if( ENABLE_SABER_DOC )
    add_subdirectory( docs )
endif()

## Global tests
ecbuild_add_test( TARGET saber_coding_norms_src
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/${PROJECT_NAME}_cpplint.py
                  ARGS --quiet --recursive ${CMAKE_SOURCE_DIR}/${PROJECT_NAME}/src )
ecbuild_add_test( TARGET saber_coding_norms_test
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/${PROJECT_NAME}_cpplint.py
                  ARGS --quiet --recursive ${CMAKE_SOURCE_DIR}/${PROJECT_NAME}/test )

## Package Config
ecbuild_install_project( NAME ${PROJECT_NAME} )

ecbuild_print_summary()
