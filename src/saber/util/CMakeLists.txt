# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set( util_src_files

tools_atlas.fypp
tools_const.F90
tools_kinds.F90
tools_repro.fypp
type_fieldset.F90
type_mpl.fypp
type_msv.fypp
type_rng.fypp
type_timer.F90

PARENT_SCOPE
)

set ( util_tmpl_files

linkedList_c.f
linkedList_i.f

PARENT_SCOPE
)
