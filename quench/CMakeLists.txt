# (C) Copyright 2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Add subdirectories
add_subdirectory( src )
add_subdirectory( mains )

if( BUILD_TESTING AND jedi-model-data_FOUND)
  # Coding norms
  ecbuild_add_test( TARGET saber_coding_norms_quench
                    TYPE SCRIPT
                    COMMAND ${CMAKE_BINARY_DIR}/bin/${PROJECT_NAME}_cpplint.py
                    ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR} )

  # Path to the saber test files
  set (jedi_model_data_tier_1_saber "${CMAKE_SOURCE_DIR}/jedi-model-data/testinput_tier_1/saber")

  # Setup QUENCH data links
  message( STATUS "Setup QUENCH data links" )
  file( REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/data )
  file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data )
  execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                   ${jedi_model_data_tier_1_saber}/landsea.nc
                   ${CMAKE_CURRENT_BINARY_DIR}/data/landsea.nc )
endif()
