# (C) Copyright 2017-2019 UCAR.
# (C) Crown Copyright 2024 Met Office
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list(APPEND util_src_files_list
BinnedFieldSetHelpers.cc
BinnedFieldSetHelpers.h
Calibration.cc
Calibration.h
HorizontalProfiles.cc
HorizontalProfiles.h
mo_cvtcoord_mod.F90
mo_netcdf_mod.F90
)

set( util_src_files ${util_src_files_list} PARENT_SCOPE )
