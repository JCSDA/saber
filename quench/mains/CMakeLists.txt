# (C) Copyright 2017-2019 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# Executables
ecbuild_add_executable( TARGET  saber_quench_error_covariance_toolbox.x
                        SOURCES quenchErrorCovarianceToolbox.cc
                        LIBS    quench
                                vader
                                ${gsibec_LIBRARIES}
                                saber )

ecbuild_add_executable( TARGET  saber_quench_process_perts.x
                        SOURCES quenchProcessPerts.cc
                        LIBS    quench
                                vader
                                ${gsibec_LIBRARIES}
                                saber )

ecbuild_add_executable( TARGET  saber_quench_convertstate.x
                        SOURCES quenchConvertState.cc
                        LIBS    quench
                                vader
                                saber )

ecbuild_add_executable( TARGET  saber_quench_coupled_error_covariance_toolbox.x
                        SOURCES quenchCoupledErrorCovarianceToolbox.cc
                        LIBS    quench
                                vader
                                ${gsibec_LIBRARIES}
                                saber )

ecbuild_add_executable( TARGET  saber_quench_subcomm_error_covariance_toolbox.x
                        SOURCES quenchSubCommErrorCovarianceToolbox.cc
                        LIBS    quench
                                vader
                                ${gsibec_LIBRARIES}
                                saber )			

# Add RPATH for Torch libraries if found
if( Torch_FOUND AND DEFINED Torch_ROOT )
    foreach(target saber_quench_error_covariance_toolbox.x saber_quench_process_perts.x saber_quench_convertstate.x saber_quench_coupled_error_covariance_toolbox.x)
        get_target_property(EXISTING_RPATH ${target} BUILD_RPATH)
        if(EXISTING_RPATH)
            set_target_properties( ${target} PROPERTIES
                BUILD_RPATH "${EXISTING_RPATH};${Torch_ROOT}/lib" )
        else()
            set_target_properties( ${target} PROPERTIES
                BUILD_RPATH "${Torch_ROOT}/lib" )
        endif()
        get_target_property(EXISTING_INSTALL_RPATH ${target} INSTALL_RPATH)
        if(EXISTING_INSTALL_RPATH)
            set_target_properties( ${target} PROPERTIES
                INSTALL_RPATH "${EXISTING_INSTALL_RPATH};${Torch_ROOT}/lib" )
        else()
            set_target_properties( ${target} PROPERTIES
                INSTALL_RPATH "${Torch_ROOT}/lib" )
        endif()
    endforeach()
endif()
