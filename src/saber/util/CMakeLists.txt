# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list(APPEND util_src_files_list
HorizontalProfiles.cc
HorizontalProfiles.h
mo_cvtcoord_mod.F90
mo_netcdf_mod.F90

# This should be removed once spack-stack uses the ATLAS version where this interface is merged
atlas_buildnode2cellconnectivity.F90
)

set( util_src_files ${util_src_files_list} PARENT_SCOPE )
