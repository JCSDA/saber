/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "saber/blocks/SaberBlockChainBase.h"
#include "saber/blocks/SaberEnsembleBlockChain.h"
#include "saber/blocks/SaberParametricBlockChain.h"

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL> void instantiateBlockChainFactory() {
  static SaberBlockChainMaker<MODEL, SaberEnsembleBlockChain>
    makerEnsembleBlockChain_("Ensemble");
  static SaberBlockChainMaker<MODEL, SaberParametricBlockChain>
    makerParametricBlockChain_("Parametric");
}

// -----------------------------------------------------------------------------

}  // namespace saber
