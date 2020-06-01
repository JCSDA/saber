# (C) Copyright 2017 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#
# SABER bundle
#

# SABER project
project( saber-bundle C CXX Fortran )

# CMake version requirement
cmake_minimum_required( VERSION 3.3.2 FATAL_ERROR )

# Define CMAKE_MODULE_PATH
set( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_CURRENT_SOURCE_DIR}/ecbuild/cmake;${CMAKE_MODULE_PATH}")
include( ecbuild_bundle )

# Initialize bundle
ecbuild_bundle_initialize()
ecbuild_requires_macro_version( 2.7 )

# Optional repositories
ecbuild_add_option( FEATURE BUILD_ECKIT
                    DEFAULT OFF
                    DESCRIPTION "download and build eckit (not needed if in a jedi container)" )

# Add repositories

# eckit
if(BUILD_ECKIT)
    ecbuild_bundle( PROJECT eckit GIT "https://github.com/JCSDA/eckit.git" TAG 1.4.0.jcsda1 )
endif()

# fckit
ecbuild_bundle( PROJECT fckit GIT "https://github.com/JCSDA/fckit.git" BRANCH develop UPDATE )

# atlas
ecbuild_bundle( PROJECT atlas GIT "https://github.com/JCSDA/atlas.git" BRANCH develop UPDATE )

# metis
ecbuild_bundle( PROJECT metis GIT "https://github.com/JCSDA/metis.git" BRANCH develop UPDATE )

# oops
ecbuild_bundle( PROJECT oops GIT "https://github.com/JCSDA/oops.git" BRANCH develop UPDATE )

# saber
ecbuild_bundle( PROJECT saber GIT "https://github.com/JCSDA/saber.git" BRANCH develop UPDATE )

# Finalize bundle
ecbuild_bundle_finalize()
