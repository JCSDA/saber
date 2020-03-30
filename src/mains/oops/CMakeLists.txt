# (C) Copyright 2017-2020 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

add_definitions(-DATLASIFIED=1)

ecbuild_add_executable( TARGET  saber_qg_estimate_parameters.x
                        SOURCES qgEstimateParams.cc
                        LIBS    saber
                      )

ecbuild_add_executable( TARGET  saber_qg_dirac.x
                        SOURCES qgDirac.cc
                        LIBS    saber
                      )

ecbuild_add_executable( TARGET  saber_qg_4dvar.x
                        SOURCES qg4DVar.cc
                        LIBS    saber
                      )
