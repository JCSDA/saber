# (C) Copyright 2017-2020 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

add_definitions(-DATLASIFIED=1)

list( APPEND interfaces_src_files
  interfacesEstimateParams.cc
)

ecbuild_add_executable( TARGET  interfaces_estimate_parameters.x
                        SOURCES ${interfaces_src_files}
                        LIBS    saber
                      )
