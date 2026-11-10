# (C) Copyright 2024 DOC/NOAA
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

if( ENABLE_MGBF )
    # MGBF Fortran library
    add_subdirectory( mgbf_lib )

    # MGBF covariance interface
    list(APPEND mgbf_src_files_list

    covariance/MGBF_Covariance.h
    covariance/MGBF_Covariance.cc
    covariance/MGBF_Covariance.interface.F90
    covariance/MGBF_Covariance.interface.h
    covariance/mgbf_covariance_mod.f90

    )

    set( mgbf_src_files

    ${mgbf_src_files_list}

    PARENT_SCOPE
    )
endif()
