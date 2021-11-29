/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_INSTANTIATESABERBLOCKFACTORY_H_
#define SABER_OOPS_INSTANTIATESABERBLOCKFACTORY_H_

#include "saber/oops/BUMP_NICAS.h"
#include "saber/oops/BUMP_PsiChiToUV.h"
#include "saber/oops/BUMP_StdDev.h"
#include "saber/oops/BUMP_VerticalBalance.h"
#include "saber/oops/GSI_RF.h"
#include "saber/oops/ID.h"
#include "saber/oops/StdDev.h"

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL>
void instantiateSaberBlockFactory() {
  // BUMP
  static SaberBlockMaker<MODEL, BUMP_NICAS<MODEL> >
              makerBUMP_NICAS_("BUMP_NICAS");
  static SaberBlockMaker<MODEL, BUMP_PsiChiToUV<MODEL> >
              makerBUMP_PsiChiToUV_("BUMP_PsiChiToUV");
  static SaberBlockMaker<MODEL, BUMP_StdDev<MODEL> >
              makerBUMP_StdDev_("BUMP_StdDev");
  static SaberBlockMaker<MODEL, BUMP_VerticalBalance<MODEL> >
              makerBUMP_VerticalBalance_("BUMP_VerticalBalance");

  // GSI_RF
  static SaberBlockMaker<MODEL, GSI_RF<MODEL> >
              makerGSI_RF_("GSI_RF");

  // Identity
  static SaberBlockMaker<MODEL, ID<MODEL> >
              makerID_("ID");

  // StdDev
  static SaberBlockMaker<MODEL, StdDev<MODEL> >
              makerStdDev_("StdDev");
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_INSTANTIATESABERBLOCKFACTORY_H_
