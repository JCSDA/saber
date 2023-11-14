# (C) Copyright 2017-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list( APPEND blocks_src_files_list
# SABER block parameters base
SaberBlockParametersBase.cc
SaberBlockParametersBase.h

# SABER central block base
SaberCentralBlockBase.cc
SaberCentralBlockBase.h

# SABER outer blocks base
SaberOuterBlockBase.cc
SaberOuterBlockBase.h

# SABER block chain base
SaberBlockChainBase.h

# SABER ensemble block chain
SaberEnsembleBlockChain.cc
SaberEnsembleBlockChain.h

# SABER outer block chain
SaberOuterBlockChain.h

# SABER parametric block chain
SaberParametricBlockChain.cc
SaberParametricBlockChain.h

# Factories
instantiateBlockChainFactory.h
)

set( blocks_src_files ${blocks_src_files_list} PARENT_SCOPE )
