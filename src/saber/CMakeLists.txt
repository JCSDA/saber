# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

add_subdirectory( bump )
add_subdirectory( external )
add_subdirectory( gaugrid )
add_subdirectory( util )

if (${OOPS_FOUND})
  add_subdirectory( oops )
endif()

PREPEND( _p_bump_src_files "bump" ${bump_src_files} )
PREPEND( _p_external_src_files "external" ${external_src_files} )
PREPEND( _p_gaugrid_src_files "gaugrid" ${gaugrid_src_files} )
PREPEND( _p_util_src_files "util" ${util_src_files} )

if (${OOPS_FOUND})
  PREPEND( _p_oops_src_files "oops" ${oops_src_files} )
endif()

list( APPEND saber_src_files
  ${_p_bump_src_files}
  ${_p_external_src_files}
  ${_p_gaugrid_src_files}
  ${_p_util_src_files}
)

if (${OOPS_FOUND})
  list( APPEND saber_src_files
    ${_p_oops_src_files}
  )
endif()

list( APPEND SABER_LIB_DEPS
  ${NETCDF_LIBRARIES} fckit eckit_mpi eckit atlas atlas_f
)

if(${METIS_FOUND})
    list( APPEND SABER_LIB_DEPS metis metis_interface )
endif()

if (${OOPS_FOUND})
    list( APPEND SABER_LIB_DEPS oops )
    if (${ENABLE_OOPS_TOYMODELS})
       list( APPEND SABER_LIB_DEPS oops qg )
    endif()
endif()

ecbuild_add_library( TARGET     saber
                     SOURCES    ${FORTRAN_BINDINGS}
                                ${saber_src_files}
                     LIBS       ${SABER_LIB_DEPS}
                     INSTALL_HEADERS LISTED
                     HEADER_DESTINATION ${INSTALL_INCLUDE_DIR}/saber
                     LINKER_LANGUAGE ${SABER_LINKER_LANGUAGE}
                   )
