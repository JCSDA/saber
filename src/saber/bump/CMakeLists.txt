# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set( bump_src_files

tools_const.F90
tools_fit.F90
tools_func.F90
tools_kinds.F90
tools_repro.F90
type_adv.F90
type_avg_blk.F90
type_avg.F90
type_bpar.F90
type_bump.F90
type_cmat_blk.F90
type_cmat.F90
type_com.F90
type_cv_blk.F90
type_cv.F90
type_diag_blk.F90
type_diag.F90
type_ens.F90
type_geom.F90
type_hdiag.F90
type_io.F90
type_lct_blk.F90
type_lct.F90
type_linop.F90
type_mesh.F90
type_minim.F90
type_mom_blk.F90
type_mom.F90
type_mpl.F90
type_msv.F90
type_nam.F90
type_nicas_blk.F90
type_nicas.F90
type_obsop.F90
type_rng.F90
type_samp.F90
type_tree.F90
type_vbal_blk.F90
type_vbal.F90

external/tools_asa007.F90
external/tools_qsort.F90
external/tools_stripack.F90

PARENT_SCOPE
)

function(generate_fortran_bindings output filename)

  set( options "" )
  set( single_value_args OUTPUT MODULE )
  set( multi_value_args "" )
  cmake_parse_arguments( _PAR "${options}" "${single_value_args}" "${multi_value_args}"  ${_FIRST_ARG} ${ARGN} )

  get_filename_component(base ${filename} NAME_WE)
  set(base_abs ${CMAKE_CURRENT_SOURCE_DIR}/${base})
  set(outfile ${CMAKE_CURRENT_BINARY_DIR}/${base}_c_binding.f90)

  if( _PAR_OUTPUT )
    set(outfile ${_PAR_OUTPUT})
  endif()
  set(${output} ${${output}} ${outfile} PARENT_SCOPE)

  if( _PAR_MODULE )
    add_custom_command(
        OUTPUT ${outfile}
        COMMAND python ${PROJECT_SOURCE_DIR}/tools/c2f.py ${CMAKE_CURRENT_SOURCE_DIR}/${filename} -o ${outfile} -m ${_PAR_MODULE}
        DEPENDS ${filename} )
  else()
    add_custom_command(
        OUTPUT ${outfile}
        COMMAND python ${PROJECT_SOURCE_DIR}/tools/c2f.py ${CMAKE_CURRENT_SOURCE_DIR}/${filename} -o ${outfile}
        DEPENDS ${filename} )
  endif()
  set_source_files_properties(${outfile} PROPERTIES GENERATED TRUE)
endfunction()
