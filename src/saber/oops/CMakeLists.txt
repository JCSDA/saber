# (C) Copyright 2017-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

if( oops_FOUND )
    list(APPEND oops_src_files_list
    # SABER block parameters base
    SaberBlockParametersBase.cc
    SaberBlockParametersBase.h

    # SABER block chain
    SaberBlockChain.cc
    SaberBlockChain.h

    # SABER central block base
    SaberCentralBlockBase.cc
    SaberCentralBlockBase.h

    # SABER outer blocks base
    SaberOuterBlockBase.cc
    SaberOuterBlockBase.h

    # ID block
    ID.h
    ID.cc

    # Ensemble block
    Ensemble.h
    Ensemble.cc

    # Hybrid block
    Hybrid.h
    Hybrid.cc

    # StdDev block
    StdDev.h
    StdDev.cc

    # Factories
    instantiateCovarFactory.h
    instantiateLocalizationFactory.h

    # Utilities
    Utilities.h
    Utilities.cc

    # Drivers
    ErrorCovariance.h
    ErrorCovarianceParameters.h
    Localization.h

    # Application
    ErrorCovarianceToolbox.h
    )
endif()

set( oops_src_files ${oops_src_files_list} PARENT_SCOPE )
