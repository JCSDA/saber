# (C) Copyright 2017-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list( APPEND oops_src_files_list
# Parameters
ErrorCovarianceParameters.h

# Drivers
ErrorCovariance.h
Localization.h

# Applications
ErrorCovarianceToolbox.h
ProcessPerts.h

# Utilities
Utilities.h
Utilities.cc

# Factories
instantiateCovarFactory.h
instantiateLocalizationFactory.h
)

set( oops_src_files ${oops_src_files_list} PARENT_SCOPE )
