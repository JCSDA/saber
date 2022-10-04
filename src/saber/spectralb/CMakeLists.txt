# (C) Crown Copyright 2022 Met Office
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list(APPEND spectralb_src_files_list
    spectralb_cvtcoord_mod.F90
    spectralb_cvtderivtypes_mod.F90
    spectralb_netcdf_mod.F90
)

if( oops_FOUND AND atlas_TRANS_FOUND )
    list(APPEND spectralb_src_files_list

    CovarianceStatistics.cc
    CovarianceStatistics.h
    CovarianceStatisticsUtils.cc
    CovarianceStatisticsUtils.h
    spectralb.h
    spectralbnointerp.h
    spectralb_covstats_interface.h
    spectralb_covstats_interface.F90
    spectralbParameters.h
    SpectralCovariance.cc
    SpectralCovariance.h
    SpectralToGauss.cc
    SpectralToGauss.h
    SpectralToGaussUV.cc
    SpectralToGaussUV.h
    SPCTRL_Cov.cc
    SPCTRL_Cov.h
    )
endif()

if( atlas_VERSION GREATER_EQUAL 0.31.0 )
    list(APPEND spectralb_src_files_list

    SpectralToGaussUV.cc
    SpectralToGaussUV.h
    )
endif()

set( spectralb_src_files

${spectralb_src_files_list}

PARENT_SCOPE
)
