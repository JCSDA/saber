# (C) Crown Copyright 2022 Met Office
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#list(vader_src_files_list
# Add future non-MO VADER blocks here
#)

# MO-specific VADER blocks
if ( oops_FOUND AND ENABLE_VADER_MO )
    list(APPEND vader_src_files_list
    AirTemperature.cc
    AirTemperature.h
    CovarianceStatisticsUtils.cc
    CovarianceStatisticsUtils.h
    DryAirDensity.cc
    DryAirDensity.h
    HydroBal.cc
    HydroBal.h
    HydrostaticExner.cc
    HydrostaticExner.h
    MoistIncrOp.cc
    MoistIncrOp.h
    MoistureControl.cc
    MoistureControl.h
    movader_covstats_interface.F90
    movader_covstats_interface.h
    )
endif()

set( vader_src_files ${vader_src_files_list} PARENT_SCOPE )
