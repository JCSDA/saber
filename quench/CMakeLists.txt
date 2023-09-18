# (C) Copyright 2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Add subdirectories
add_subdirectory( src )
add_subdirectory( mains )

# Coding norms
ecbuild_add_test( TARGET saber_coding_norms_quench
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/${PROJECT_NAME}_cpplint.py
                  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR} )

# Setup QUENCH data links
message( STATUS "Setup QUENCH data links" )
file( REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/data )
file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data )
execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                 ${CMAKE_CURRENT_SOURCE_DIR}/data/landsea.nc
                 ${CMAKE_CURRENT_BINARY_DIR}/data/landsea.nc )
