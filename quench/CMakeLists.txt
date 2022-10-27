# (C) Copyright 2022 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

add_subdirectory( src )
add_subdirectory( mains )

ecbuild_add_test( TARGET saber_coding_norms_quench
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/oops_cpplint.py
                  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR} )
