# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list(APPEND interpolation_src_files_list
# ATLAS interpolation wrapper
AtlasInterpWrapper.cc
AtlasInterpWrapper.h

# Rescaling block for Atlas interpolation wrapper
Rescaling.cc
Rescaling.h

# Interpolation block
Geometry.cc
Geometry.h
Interpolation.cc
Interpolation.h

# VertProj block
VertProj.cc
VertProj.h
)

if( atlas_TRANS_FOUND OR atlas_ECTRANS_FOUND )
    list(APPEND interpolation_src_files_list
        # ATLAS gauss to cubed-sphere-dual saber block
        GaussToCS.cc
        GaussToCS.h
    )
endif()

set( interpolation_src_files

${interpolation_src_files_list}

PARENT_SCOPE
)
