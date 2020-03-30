# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# SABER
################################################################################

cmake_minimum_required( VERSION 3.3.2 FATAL_ERROR )

project( saber C CXX Fortran )

set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH})

set( CMAKE_DIRECTORY_LABELS "saber" )

set( ECBUILD_DEFAULT_BUILD_TYPE Release )
set( ENABLE_OS_TESTS           OFF CACHE BOOL "Disable OS tests" FORCE )
set( ENABLE_LARGE_FILE_SUPPORT OFF CACHE BOOL "Disable testing of large file support" FORCE )
set( BUILD_BUMP_DOC OFF )

include( ecbuild_system NO_POLICY_SCOPE )

ecbuild_requires_macro_version( 2.5 )

################################################################################
# Project
################################################################################

ecbuild_declare_project()

ecbuild_enable_fortran( REQUIRED )
ecbuild_add_cxx11_flags()

set( SABER_LINKER_LANGUAGE CXX )

################################################################################
# Standard compiler flags
################################################################################

include( saber_compiler_flags )
include( saber_extra_macros )

################################################################################
# Dependencies
################################################################################

# Boost
include_directories( ${Boost_INCLUDE_DIR} )

# OpenMP
ecbuild_add_option( FEATURE OMP
                    DESCRIPTION "Support for OpenMP threaded parallelism"
                    REQUIRED_PACKAGES "OMP COMPONENTS CXX Fortran" )
if( HAVE_OMP )
  ecbuild_enable_omp()
else()
  ecbuild_enable_ompstubs()
endif()

# NetCDF
find_package( NetCDF REQUIRED COMPONENTS Fortran )
include_directories( ${NETCDF_INCLUDE_DIRS} )

# eckit and MPI
ecbuild_use_package( PROJECT eckit VERSION 1.1.0 REQUIRED )
include_directories( ${ECKIT_INCLUDE_DIRS} )
if ( NOT ECKIT_HAVE_MPI )
  ecbuild_error("SABER requested to build with MPI but eckit wasn't built with MPI: ECKIT_HAVE_MPI = ${ECKIT_HAVE_MPI}")
endif()

# fckit
ecbuild_use_package( PROJECT fckit VERSION 0.6.3 REQUIRED )
include_directories( ${FCKIT_INCLUDE_DIRS} )

# atlas
ecbuild_use_package( PROJECT atlas VERSION 0.18.1 REQUIRED )
include_directories( ${ATLAS_INCLUDE_DIRS} )

# git
find_package( Git )

# oops
ecbuild_use_package( PROJECT oops VERSION 0.2.1 REQUIRED )
include_directories( ${OOPS_INCLUDE_DIRS} )

################################################################################
# Export package info
################################################################################

list( APPEND SABER_TPLS eckit eckit_mpi fckit atlas atlas_f oops qg )

set( SABER_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src
                        ${CMAKE_CURRENT_BINARY_DIR}/src
                        ${CMAKE_Fortran_MODULE_DIRECTORY} )

set( SABER_LIBRARIES saber )

get_directory_property( SABER_DEFINITIONS COMPILE_DEFINITIONS )

foreach( _tpl ${SABER_TPLS} )
  string( TOUPPER ${_tpl} TPL )
  list( APPEND SABER_EXTRA_DEFINITIONS   ${${TPL}_DEFINITIONS}  ${${TPL}_TPL_DEFINITIONS}  )
  list( APPEND SABER_EXTRA_INCLUDE_DIRS  ${${TPL}_INCLUDE_DIRS} ${${TPL}_TPL_INCLUDE_DIRS} )
  list( APPEND SABER_EXTRA_LIBRARIES     ${${TPL}_LIBRARIES}    ${${TPL}_TPL_LIBRARIES}    )
endforeach()

################################################################################
# Sources
################################################################################

include_directories( ${SABER_INCLUDE_DIRS}
                     ${SABER_EXTRA_INCLUDE_DIRS}
                     ${CMAKE_CURRENT_SOURCE_DIR}/src )

add_subdirectory( doc )
add_subdirectory( src )
add_subdirectory( test )
add_subdirectory( tools )

if(ECBUILD_INSTALL_FORTRAN_MODULES)
  install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/${CMAKE_CFG_INTDIR} DESTINATION ${INSTALL_INCLUDE_DIR} )
endif()

################################################################################
# Finalise configuration
################################################################################

# prepares a tar.gz of the sources and/or binaries
ecbuild_install_project( NAME saber )

# print the summary of the configuration
ecbuild_print_summary()

ecbuild_add_test( TARGET ufo_coding_norms
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/cpplint.py
                  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR}
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
