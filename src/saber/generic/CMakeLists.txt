# (C) Copyright 2017-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list(APPEND generic_src_files_list
# Duplicate Variables
DuplicateVariables.cc
DuplicateVariables.h

# ID block
ID.h
ID.cc

# Orographic interpolation block
OrographicInterp.h
OrographicInterp.cc

# Residual fields block
ResidualFields.h
ResidualFields.cc

# Shadow levels block
ShadowLevels.h
ShadowLevels.cc

# StdDev block
StdDev.h
StdDev.cc

# Vertical localization
VertLoc.cc
VertLoc.h

# Interpolated vertical localization
VertLocInterp.cc
VertLocInterp.h
)

# Write selected fields to a netCDF file (depends of ECTRANS)
if( atlas_TRANS_FOUND OR atlas_ECTRANS_FOUND )
  list(APPEND generic_src_files_list
    WriteFields.h
    WriteFields.cc
  )
endif()

set( generic_src_files ${generic_src_files_list} PARENT_SCOPE )
