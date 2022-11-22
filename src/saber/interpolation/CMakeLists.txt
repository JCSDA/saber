# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list(APPEND interpolation_src_files_list
# BUMP interpolator
interpolatorbump_mod.fypp
)

if( oops_FOUND )
    list(APPEND interpolation_src_files_list
    # ATLAS interpolation wrapper
    AtlasInterpWrapper.h
    AtlasInterpWrapper.cc
    )

    if( atlas_TRANS_FOUND )
        list(APPEND interpolation_src_files_list
            # ATLAS gauss to cubed-sphere-dual saber block
            GaussToCS.h
            GaussToCS.cc
        )
    endif()
endif()


set( interpolation_src_files

${interpolation_src_files_list}

PARENT_SCOPE
)
