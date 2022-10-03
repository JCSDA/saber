# (C) Crown Copyright 2022 Met Office
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list(APPEND spectralb_src_files_list
    spectralb_cvtcoord_mod.F90
    spectralb_cvtderivtypes_mod.F90
    spectralb_netcdf_mod.F90
)

if( atlas_TRANS_FOUND )
    list(APPEND spectralb_src_files_list

    CovarianceStatistics.cc
    CovarianceStatistics.h
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
    SPCTRL_Cov.cc
    SPCTRL_Cov.h
    SPTOUV.cc
    SPTOUV.h
    )
endif()

set( spectralb_src_files

${spectralb_src_files_list}

PARENT_SCOPE
)
