# (C) Crown Copyright 2022 Met Office
# (C) Copyright 2024 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list(APPEND vader_src_files_list
# MO-specific VADER blocks
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
VarianceAccumulationUtils.cc
VarianceAccumulationUtils.h
WriteVariances.cc
WriteVariances.h
# Generic VADER blocks
DefaultCookbook.h
VaderBlock.cc
VaderBlock.h
)

set( vader_src_files ${vader_src_files_list} PARENT_SCOPE )
