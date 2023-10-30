# (C) Copyright 2017-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list(APPEND generic_src_files_list
# ID block
ID.h
ID.cc

# Duplicate Variables
DuplicateVariables.cc
DuplicateVariables.h

# Ensemble block
Ensemble.h
Ensemble.cc

# Hybrid block
Hybrid.h
Hybrid.cc

# StdDev block
StdDev.h
StdDev.cc

# Write selected fields to a netCDF file
WriteFields.h
WriteFields.cc
)

set( generic_src_files ${generic_src_files_list} PARENT_SCOPE )
