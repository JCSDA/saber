# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set( util_src_files

tools_const.F90
tools_kinds.F90
tools_repro.F90
type_mpl.F90
type_msv.F90
type_rng.F90

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
