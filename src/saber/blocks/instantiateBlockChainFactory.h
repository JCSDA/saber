/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "saber/blocks/BlockChainBase.h"
#include "saber/blocks/EnsembleBlockChain.h"
#include "saber/blocks/HybridBlockChain.h"
#include "saber/blocks/ParametricBlockChain.h"
#if defined(GSIBEC_FOUND)
#include "saber/gsi/GSIBlockChain.h"
#endif

namespace saber {

// -----------------------------------------------------------------------------
template <typename MODEL>
void instantiateBlockChainFactory() {
  static BlockChainMaker<MODEL, EnsembleBlockChain>
    makerEnsembleBlockChain_("ensemble");
  static BlockChainMaker<MODEL, HybridBlockChain<MODEL>>
    makerHybridBlockChain_("hybrid");
  static BlockChainMaker<MODEL, ParametricBlockChain>
    makerParametricBlockChain_("parametric");
#if defined(GSIBEC_FOUND)
  static BlockChainMaker<MODEL, gsi::GSIBlockChain>
    makerGSIBlockChain_("gsi hybrid covariance");
#endif
}

// -----------------------------------------------------------------------------

}  // namespace saber
