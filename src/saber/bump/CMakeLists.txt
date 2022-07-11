# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list(APPEND bump_src_files_list

tools_fit.fypp
tools_func.c
tools_func.fypp
tools_gc99.fypp
type_avg.fypp
type_avg_blk.fypp
type_bpar.fypp
type_bump.fypp
type_bump.h
type_bump_interface.F90
type_cmat.fypp
type_cmat_blk.fypp
type_com.fypp
type_cv.fypp
type_cv_blk.fypp
type_cv_cmp.fypp
type_diag.fypp
type_diag_blk.fypp
type_ens.fypp
type_geom.fypp
type_hdiag.fypp
type_io.fypp
type_linop.fypp
type_mesh.fypp
type_minim.fypp
type_mom.fypp
type_mom_blk.fypp
type_nam.fypp
type_nicas.fypp
type_nicas_blk.fypp
type_nicas_cmp.fypp
type_samp.fypp
type_tree.fypp
type_var.fypp
type_vbal.fypp
type_vbal_blk.fypp
type_wind.fypp

)

if( oops_FOUND )
    list(APPEND bump_src_files_list

    BUMP.h
    BUMP_NICAS.h
    BUMP_PsiChiToUV.h
    BUMP_StdDev.h
    BUMP_VerticalBalance.h

    )
endif()

set( bump_src_files

${bump_src_files_list}

PARENT_SCOPE
)
