# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list(APPEND util_src_files_list

tools_atlas.fypp
tools_const.fypp
tools_kinds.fypp
tools_netcdf.fypp
tools_repro.fypp
type_fieldset.fypp
type_mpl.fypp
type_msv.fypp
type_rng.fypp

)

if( ENABLE_SABER_INSTRUMENTATION )
    list(APPEND util_src_files_list

    type_probe.fypp

    )
endif()

if ( NOT oops_FOUND )

    list(APPEND util_src_files_list

    logger_mod.fypp

    )

endif()

set( util_src_files

${util_src_files_list}

PARENT_SCOPE
)
