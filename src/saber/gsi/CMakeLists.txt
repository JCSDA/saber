# (C) Copyright 2022 United States Government as represented by the Administrator of the National
#     Aeronautics and Space Administration
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

if( gsibec_FOUND )
    list(APPEND gsi_src_files_list
    # Covariance block
    covariance/Covariance.cc
    covariance/Covariance.h
    covariance/Covariance.interface.F90
    covariance/Covariance.interface.h
    covariance/gsi_covariance_mod.f90

    # Grid
    grid/Grid.cc
    grid/Grid.h
    grid/Grid.interface.F90
    grid/Grid.interface.h
    grid/gsi_grid_mod.f90

    # Interpolation block
    interpolation/Interpolation.cc
    interpolation/Interpolation.h

    # Unstructured interpolation code ported from oops (until new interp code can be used)
    interpolation/unstructured_interp/saber_unstructured_interpolation_mod.F90
    interpolation/unstructured_interp/UnstructuredInterpolation.cc
    interpolation/unstructured_interp/UnstructuredInterpolation.h
    interpolation/unstructured_interp/UnstructuredInterpolation.interface.F90
    interpolation/unstructured_interp/UnstructuredInterpolation.interface.h

    # Utilities
    utils/GSIParameters.h
    utils/gsi_utils_mod.f90
    )
endif()

set( gsi_src_files ${gsi_src_files_list} PARENT_SCOPE )
