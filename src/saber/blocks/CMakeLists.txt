# (C) Copyright 2017-2021 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

list( APPEND blocks_src_files_list
# SABER block parameters base
BlockParametersBase.cc
BlockParametersBase.h

# SABER central block base
CentralBlockBase.cc
CentralBlockBase.h

# SABER outer blocks base
OuterBlockBase.cc
OuterBlockBase.h

# SABER block chain base
BlockChainBase.h

# SABER central block wrapper
CentralBlockWrapper.cc
CentralBlockWrapper.h

# SABER ensemble block chain
EnsembleBlockChain.cc
EnsembleBlockChain.h

# SABER hybrid block chain
HybridBlockChain.h

# SABER outer block chain
OuterBlockChain.cc
OuterBlockChain.h

# SABER parametric block chain
ParametricBlockChain.cc
ParametricBlockChain.h

# Factories
instantiateBlockChainFactory.h
)

set( blocks_src_files ${blocks_src_files_list} PARENT_SCOPE )
