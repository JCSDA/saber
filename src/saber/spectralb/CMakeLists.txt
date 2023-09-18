# (C) Crown Copyright 2022 Met Office
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

if( atlas_TRANS_FOUND OR atlas_ECTRANS_FOUND )
    list(APPEND spectralb_src_files_list
    CovarianceStatistics.cc
    CovarianceStatistics.h
    CovarianceStatisticsUtils.cc
    CovarianceStatisticsUtils.h
    SpectralAnalyticalFilter.cc
    SpectralAnalyticalFilter.h
    spectralb_covstats_interface.h
    spectralb_covstats_interface.F90
    spectralbParameters.h
    MOSpectralCovariance.cc
    MOSpectralCovariance.h
    SpectralCovariance.cc
    SpectralCovariance.h
    SpectralToGauss.cc
    SpectralToGauss.h
    SpectralToSpectral.cc
    SpectralToSpectral.h
    SqrtOfSpectralCovariance.cc
    SqrtOfSpectralCovariance.h
    )

    if( vader_FOUND )
          list(APPEND spectralb_src_files_list
          GaussUVToGP.cc
          GaussUVToGP.h
          HydrostaticPressure.cc
          HydrostaticPressure.h
          )
    endif()
endif()

set( spectralb_src_files ${spectralb_src_files_list} PARENT_SCOPE )
