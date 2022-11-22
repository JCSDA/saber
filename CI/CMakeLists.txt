# (C) Copyright 2017 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#
# SABER bundle
#

cmake_minimum_required( VERSION 3.12 FATAL_ERROR )

project( saber-bundle VERSION 1.0.0 LANGUAGES CXX Fortran )

set( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_CURRENT_SOURCE_DIR}/ecbuild/cmake;${CMAKE_MODULE_PATH}")
include( ecbuild_bundle )

set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )

ecbuild_bundle_initialize()

ecbuild_requires_macro_version( 2.7 )

ecbuild_add_option( FEATURE OMP
                    DEFAULT ON
                    DESCRIPTION "Use OpenMP" )

#ecbuild_bundle( PROJECT jedicmake GIT "https://github.com/JCSDA-internal/jedi-cmake.git" )
#include( jedicmake/cmake/Functions/git_functions.cmake  )
include( /opt/view/share/jedicmake/Functions/git_functions.cmake )

ecbuild_bundle( PROJECT oops GIT "https://github.com/JCSDA-internal/oops.git" )
ecbuild_bundle( PROJECT vader GIT "https://github.com/JCSDA-internal/vader.git" )
ecbuild_bundle( PROJECT gsibec GIT "https://github.com/GEOS-ESM/GSIbec.git" )
ecbuild_bundle( PROJECT saber GIT "https://github.com/JCSDA-internal/saber.git" )

ecbuild_bundle_finalize()

include(cmake/cdash-integration.cmake)
include(CTest)
