# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list(APPEND interpolation_src_files_list
interpolatorbump_mod.fypp
)

if( oops_FOUND )
    list(APPEND interpolation_src_files_list
    AtlasInterpWrapper.h
    )
endif()

set( interpolation_src_files ${interpolation_src_files_list} PARENT_SCOPE )
