# (C) Crown Copyright 2022 Met Office
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#list(vader_src_files_list
# Add future non-MO VADER blocks here
#)

# MO-specific VADER blocks
if ( oops_FOUND AND vader_FOUND )
    list(APPEND vader_src_files_list
    AirTemperature.cc
    AirTemperature.h
    CovarianceStatisticsUtils.cc
    CovarianceStatisticsUtils.h
    DryAirDensity.cc
    DryAirDensity.h
    GpToHp.cc
    GpToHp.h
    HpToHexner.cc
    HpToHexner.h
    HpHexnerToPExner.cc
    HpHexnerToPExner.h
    HydroBal.cc
    HydroBal.h
    HydrostaticExner.cc
    HydrostaticExner.h
    MoistIncrOp.cc
    MoistIncrOp.h
    MoistIncrOpParameters.h
    MoistureControl.cc
    MoistureControl.h
    movader_covstats_interface.F90
    movader_covstats_interface.h
    SuperMoistIncrOp.cc
    SuperMoistIncrOp.h
    PressureParameters.h
    VertLoc.cc
    VertLoc.h
    VertLocInterp.cc
    VertLocInterp.h
    )
endif()

if ( vader_FOUND )
    list(APPEND vader_src_files_list

    VaderBlock.cc
    VaderBlock.h
    )
endif()

set( vader_src_files ${vader_src_files_list} PARENT_SCOPE )
