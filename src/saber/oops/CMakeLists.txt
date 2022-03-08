# (C) Copyright 2017-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set( oops_src_files

# SABER blocks base
SaberBlockBase.h
SaberBlockParametersBase.h

# BUMP blocks
BUMP.h
BUMP_NICAS.h
BUMP_PsiChiToUV.h
BUMP_StdDev.h
BUMP_VerticalBalance.h

# GSI block
GSI.h
GSI_RF.h

# ID block
ID.h

# StdDev block
StdDev.h

# Factories
instantiateCovarFactory.h
instantiateLocalizationFactory.h
instantiateSaberBlockFactory.h

# Drivers
ErrorCovariance.h
ErrorCovarianceTraining.h
Localization.h
SaberBlockTest.h

PARENT_SCOPE
)
