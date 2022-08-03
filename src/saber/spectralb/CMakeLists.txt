# (C) Crown Copyright 2022 Met Office
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

if( atlas_TRANS_FOUND )
    list(APPEND spectralb_src_files_list

    CovarianceStatistics.h
    CovarianceStatisticsUtils.h
    spectralb.h
    spectralb_covstats_interface.h
    spectralb_covstats_interface.F90
    spectralb_cvtcoord_mod.F90
    spectralb_cvtderivtypes_mod.F90
    spectralb_netcdf_mod.F90
    spectralbParameters.h
    SPCTRL_Cov.h

    )
endif()

set( spectralb_src_files

${spectralb_src_files_list}

PARENT_SCOPE
)
