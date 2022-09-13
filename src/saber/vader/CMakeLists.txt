# (C) Crown Copyright 2022 Met Office
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#list(vader_src_files_list

# Add future non-MO VADER blocks here

#)

# MO-specific VADER blocks
if ( ENABLE_VADER_MO )
    list(APPEND vader_src_files_list

    AirTemperatureSaberBlock.cc
    AirTemperatureSaberBlock.h
    CovarianceStatisticsUtils.cc
    CovarianceStatisticsUtils.h
    DryAirDensitySaberBlock.cc
    DryAirDensitySaberBlock.h
    HydroBalSaberBlock.cc
    HydroBalSaberBlock.h
    HydrostaticExnerParameters.h
    HydrostaticExnerSaberBlock.cc
    HydrostaticExnerSaberBlock.h
    MoistIncrOpSaberBlock.cc
    MoistIncrOpSaberBlock.h
    MoistureControlParameters.h
    MoistureControlSaberBlock.cc
    MoistureControlSaberBlock.h
    movader_covstats_interface.F90
    movader_covstats_interface.h
    )
endif()

set( vader_src_files

${vader_src_files_list}

PARENT_SCOPE
)
