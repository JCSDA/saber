# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list(APPEND bump_lib_src_files_list
BUMP.cc
BUMP.h
interpolatorbump_mod.fypp
tools_asa007.fypp
tools_atlas.fypp
tools_const.fypp
tools_fit.fypp
tools_func.c
tools_func.fypp
tools_gc99.fypp
tools_kinds.fypp
tools_log.cc
tools_log.fypp
tools_netcdf.fypp
tools_qsort.fypp
tools_repro.fypp
tools_ssrfpack.fypp
tools_stripack.fypp
tools_wrfda.fypp
type_avg.fypp
type_avg_blk.fypp
type_bump.h
type_bump.fypp
type_bump_interface.F90
type_bump_parameters.cc
type_bump_parameters.h
type_cmat.fypp
type_cmat_blk.fypp
type_com.fypp
type_cv.fypp
type_cv_blk.fypp
type_cv_cmp.fypp
type_diag.fypp
type_diag_blk.fypp
type_ens.fypp
type_fieldset.fypp
type_geom.fypp
type_gsi.fypp
type_hdiag.fypp
type_io.fypp
type_linop.fypp
type_mesh.fypp
type_minim.fypp
type_mom.fypp
type_mom_blk.fypp
type_mpl.fypp
type_msv.fypp
type_nam.fypp
type_nicas.fypp
type_nicas_blk.fypp
type_nicas_cmp.fypp
type_rng.fypp
type_samp.fypp
type_tree.fypp
type_var.fypp
type_vbal.fypp
type_vbal_blk.fypp
type_wind.fypp
)

if( ENABLE_SABER_INSTRUMENTATION )
    list(APPEND bump_lib_src_files_list
    type_probe.fypp
    )
endif()

set( bump_lib_src_files ${bump_lib_src_files_list} PARENT_SCOPE )
