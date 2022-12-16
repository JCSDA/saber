# (C) Crown Copyright 2022 Met Office
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

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
    SPCTRL_Cov.cc
    SPCTRL_Cov.h
    )

    if( vader_FOUND AND atlas_VERSION GREATER_EQUAL 0.31.0 )
          list(APPEND spectralb_src_files_list
          GaussUVToGP.cc
          GaussUVToGP.h
          )
    endif()
endif()

set( spectralb_src_files ${spectralb_src_files_list} PARENT_SCOPE )
