# (C) Copyright 2017-2019 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 
# VADER libraries
if( vader_FOUND )
    set( vader_LIBRARIES "vader" )
else()
    set( vader_LIBRARIES "" )
endif()

# Executables
ecbuild_add_executable( TARGET  saber_quench_randomization.x
                        SOURCES quenchRandomization.cc
                        LIBS    quench
                                ${vader_LIBRARIES}
                                ${gsibec_LIBRARIES}
                                saber )

ecbuild_add_executable( TARGET  saber_quench_error_covariance_training.x
                        SOURCES quenchErrorCovarianceTraining.cc
                        LIBS    quench
                                ${vader_LIBRARIES}
                                ${gsibec_LIBRARIES}
                                saber )

ecbuild_add_executable( TARGET  saber_quench_saber_block_test.x
                        SOURCES quenchSaberBlockTest.cc
                        LIBS    quench
                                ${vader_LIBRARIES}
                                ${gsibec_LIBRARIES}
                                saber )

ecbuild_add_executable( TARGET  saber_quench_dirac.x
                        SOURCES quenchDirac.cc
                        LIBS    quench
                                ${vader_LIBRARIES}
                                ${gsibec_LIBRARIES}
                                saber )

ecbuild_add_executable( TARGET  saber_quench_convertstate.x
                        SOURCES quenchConvertState.cc
                        LIBS    quench
                                ${vader_LIBRARIES}
                                saber )
