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

## RPATH settings - ensure executables can find shared libraries during build and install
set( CMAKE_SKIP_BUILD_RPATH FALSE )
set( CMAKE_BUILD_WITH_INSTALL_RPATH FALSE )
set( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )
# Add build lib directory to RPATH
set( CMAKE_BUILD_RPATH "${CMAKE_BINARY_DIR}/lib" )
list( APPEND CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" )

## Configuration options
option( ENABLE_OFFLINE_CODECOV "Enable SABER offline code coverage" OFF )
option( ENABLE_SABER_DOC "Enable SABER documentation" OFF )
option( ENABLE_MKL "Use MKL for LAPACK implementation (if available)" ON )
option( ENABLE_BUMP "Compile BUMP and run BUMP-based tests" ON )
option( ENABLE_BUMP_INSTRUMENTATION "Enable BUMP instrumentation" OFF )
option( ENABLE_QUENCH "Compile QUENCH and run QUENCH-based tests" ON )
include( ${PROJECT_NAME}_compiler_flags )

## Dependencies
# Required
if( ENABLE_QUENCH )
    find_package( jedicmake REQUIRED )  # Prefer find modules from jedi-cmake
endif()
if(OPENMP)
    find_package( OpenMP REQUIRED COMPONENTS CXX Fortran )
endif()
find_package( MPI REQUIRED COMPONENTS Fortran )
find_package( NetCDF REQUIRED COMPONENTS C Fortran )
find_package( eckit 1.17.1 REQUIRED COMPONENTS MPI )
find_package( fckit 0.9.5 REQUIRED )
if(OpenMP_FOUND)
    find_package( atlas 0.35.0 REQUIRED COMPONENTS OMP_Fortran TESSELATION)
else()
    find_package( atlas 0.35.0 REQUIRED COMPONENTS TESSELATION)
endif()
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

# Flag to signify Atlas supports spherical vector interpolation.
if( atlas_VERSION VERSION_GREATER_EQUAL "0.44.0" )
    add_definitions(-DATLAS_SUPPORTS_SPHERICAL_VECTOR_INTERP_META)
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
find_package( ECTRANS 1.6.0 COMPONENTS transi double QUIET ) # TODO(Benjamin): update tag once the ectrans PR is merged and a new tag is created
find_package( FFTW 3.3.8 QUIET )
find_package( ip 5 QUIET )
find_package( gsibec 1.4.2 QUIET )

# Torch
set(Torch_FOUND FALSE CACHE BOOL "Torch not found")
if(DEFINED Torch_ROOT AND NOT Torch_ROOT STREQUAL "")
  set(Torch_FOUND TRUE CACHE BOOL "Torch found" FORCE)
  # Temporary fix until Torch is properly installed as a cmake package
  include_directories("${Torch_ROOT}/include")
  include_directories("${Torch_ROOT}/include/torch/csrc/api/include")
  # Use full paths to libraries instead of link_directories
  set(TORCH_LIBRARIES
      "${Torch_ROOT}/lib/libtorch.so"
      "${Torch_ROOT}/lib/libtorch_cpu.so"
      "${Torch_ROOT}/lib/libc10.so")
endif()

# Optional SABER blocks
if( FFTW_FOUND OR ( ECTRANS_FOUND AND ectrans_HAVE_ETRANS AND ectrans_HAVE_TRANSI ) )
    message( STATUS "SABER block Bifourier is enabled" )
else()
    message( STATUS "SABER block Bifourier is NOT enabled" )
endif()
if( FFTW_FOUND )
    message( STATUS "SABER block FastLAM spectral layer is enabled" )
else()
    message( STATUS "SABER block FastLAM spectral layer is NOT enabled" )
endif()
if( gsibec_FOUND )
    message( STATUS "SABER block GSI is enabled" )
else()
    message( STATUS "SABER block GSI is NOT enabled" )
endif()
if( atlas_TRANS_FOUND OR atlas_ECTRANS_FOUND )
    message( STATUS "SABER block SPECTRALB is enabled" )
else()
    message( STATUS "SABER block SPECTRALB is NOT enabled" )
endif()

## Sources
add_subdirectory( src/saber )
add_subdirectory( tools )

## BUMP instrumentation
if( ENABLE_BUMP )
    if( ENABLE_BUMP_INSTRUMENTATION )
        message( STATUS "BUMP instrumentation is enabled" )
        add_definitions(-DENABLE_BUMP_INSTRUMENTATION=1)
    else()
        message( STATUS "BUMP instrumentation is not enabled" )
        add_definitions(-DENABLE_BUMP_INSTRUMENTATION=0)
    endif()
endif()

if( ENABLE_SABER_DOC )
    message( STATUS "SABER documentation is enabled" )
    add_subdirectory( docs )
else()
    message( STATUS "SABER documentation is not enabled" )
endif()

if( BUILD_TESTING )
  # Test data
  find_package( jedi-model-data 1.0.0 QUIET )

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

  # If test data repo found build the tests
  if( jedi-model-data_FOUND )
      message( STATUS "jedi-model-data found; building saber tests" )
      add_subdirectory( test )
  else()
      message( STATUS "jedi-model-data not found; skipping saber tests" )
  endif()

  ## Global tests
  ecbuild_add_test( TARGET saber_coding_norms_src
                    TYPE SCRIPT
                    COMMAND ${CMAKE_BINARY_DIR}/bin/${PROJECT_NAME}_cpplint.py
                    ARGS --quiet --recursive ${CMAKE_SOURCE_DIR}/${PROJECT_NAME}/src )
endif()
if( ENABLE_QUENCH )
  add_subdirectory( quench )
endif()

## Package Config
ecbuild_install_project( NAME ${PROJECT_NAME} )

ecbuild_print_summary()
