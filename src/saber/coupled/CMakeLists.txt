# (C) Copyright 2025- UCAR.
# (C) Copyright 2025- NOAA/EMC.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list( APPEND coupled_src_files_list
instantiateCoupledCovarFactory.h

CoupledErrorCovariance.h

CoupledErrorCovarianceParameters.h
)

set( coupled_src_files ${coupled_src_files_list} PARENT_SCOPE )
