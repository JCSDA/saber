# (C) Copyright 2017-2020 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# SABER
################################################################################

cmake_minimum_required( VERSION 3.12 )

project( saber VERSION 1.10.0 LANGUAGES CXX Fortran )

option(OPENMP "Build saber with OpenMP support" ON)

## Ecbuild integration
find_package( ecbuild QUIET )
include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_declare_project()
list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
set( CMAKE_DIRECTORY_LABELS ${PROJECT_NAME} )

## Configuration options
option( ENABLE_SABER_INSTRUMENTATION "Enable SABER instrumentation" OFF )
option( ENABLE_OFFLINE_CODECOV "Enable SABER offline code coverage" OFF )
option( ENABLE_SABER_DOC "Enable SABER documentation" OFF )
option( ENABLE_MKL "Use MKL for LAPACK implementation (if available)" ON )
option( ENABLE_QUENCH "Compile QUENCH and run QUENCH-based tests" ON )
include( ${PROJECT_NAME}_compiler_flags )

## Dependencies
# Required
find_package( jedicmake REQUIRED )  # Prefer find modules from jedi-cmake
if(OPENMP)
  find_package( OpenMP REQUIRED COMPONENTS CXX Fortran )
endif()
find_package( MPI REQUIRED COMPONENTS Fortran )
find_package( NetCDF REQUIRED COMPONENTS C Fortran )
find_package( eckit 1.17.1 REQUIRED COMPONENTS MPI )
find_package( fckit 0.9.5 REQUIRED )
find_package( atlas 0.35.0 REQUIRED COMPONENTS OMP_Fortran TESSELATION)
# TODO(Benjamin): check when ATLAS PR #215 is merged into a tagged version to update the version number
if( atlas_VERSION VERSION_GREATER "0.38.1" )
    set( ATLAS_REGIONAL_INTERP ON )
else()
    set( ATLAS_REGIONAL_INTERP OFF )
endif()
if( DEFINED ENV{ATLAS_REGIONAL_INTERP} )
    set( ATLAS_REGIONAL_INTERP $ENV{ATLAS_REGIONAL_INTERP} )
endif()

if( atlas_VERSION VERSION_GREATER "0.41" )
    add_definitions(-DATLAS_MAKE_SPARSE)
endif()

if( ENABLE_MKL )
    find_package( MKL )
endif()
if( ENABLE_MKL AND MKL_FOUND )
    set( LAPACK_LIBRARIES ${MKL_LIBRARIES} )
else()
    find_package( LAPACK REQUIRED )
endif()

# Required
find_package( oops 1.10.0 REQUIRED )
find_package( vader 1.7.0 REQUIRED )

# Optionals
find_package( FFTW 3.3.8 QUIET )
find_package( gsibec 1.2.1 QUIET )
find_package( eccodes 2.24 QUIET )
if( eccodes_FOUND )
    set( ECCODES_LIBRARIES eccodes )
    add_definitions(-DECCODES_FOUND=1)
endif()

# Optional SABER blocks
if( gsibec_FOUND )
    find_package( sp QUIET )
    message( STATUS "SABER block GSI is enabled" )
else()
    message( STATUS "SABER block GSI is NOT enabled" )
endif()
if( FFTW_FOUND )
    message( STATUS "SABER block FastLAM spectral layer is enabled" )
else()
    message( STATUS "SABER block FastLAM spectral layer is NOT enabled" )
endif()
if( atlas_TRANS_FOUND OR atlas_ECTRANS_FOUND )
    message( STATUS "SABER block SPECTRALB is enabled" )
else()
    message( STATUS "SABER block SPECTRALB is NOT enabled" )
endif()

## SABER instrumentation
if( ENABLE_SABER_INSTRUMENTATION )
    message( STATUS "SABER instrumentation is enabled" )
    add_definitions(-DENABLE_SABER_INSTRUMENTATION=1)
else()
    message( STATUS "SABER instrumentation is not enabled" )
    add_definitions(-DENABLE_SABER_INSTRUMENTATION=0)
endif()

## Code coverage
if ( ENABLE_OFFLINE_CODECOV )
    message( STATUS "Offline code coverage is enabled" )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    set( CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -fprofile-arcs -ftest-coverage")
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fprofile-arcs -ftest-coverage")
    set( CMAKE_Fortran_LINK_FLAGS "${CMAKE_Fortran_LINK_FLAGS} -fprofile-arcs -ftest-coverage")
else()
    message( STATUS "Offline code coverage is not enabled" )
endif()

## Sources
if( ENABLE_QUENCH )
  add_subdirectory( quench )
endif()
add_subdirectory( src/saber )
add_subdirectory( tools )
add_subdirectory( test )
if( ENABLE_SABER_DOC )
    message( STATUS "SABER documentation is enabled" )
    add_subdirectory( docs )
else()
    message( STATUS "SABER documentation is not enabled" )
endif()

## Global tests
ecbuild_add_test( TARGET saber_coding_norms_src
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/${PROJECT_NAME}_cpplint.py
                  ARGS --quiet --recursive ${CMAKE_SOURCE_DIR}/${PROJECT_NAME}/src )

## Package Config
ecbuild_install_project( NAME ${PROJECT_NAME} )

ecbuild_print_summary()
