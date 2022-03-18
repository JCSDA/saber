/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_INSTANTIATESABERBLOCKFACTORY_H_
#define SABER_OOPS_INSTANTIATESABERBLOCKFACTORY_H_

#include "saber/bump/BUMP_NICAS.h"
#include "saber/bump/BUMP_PsiChiToUV.h"
#include "saber/bump/BUMP_StdDev.h"
#include "saber/bump/BUMP_VerticalBalance.h"
#include "saber/gsi/covariance/GSI_Covariance.h"
#include "saber/gsi/interpolation/GSI_Interpolation.h"
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

  // GSI operators
  static SaberBlockMaker<MODEL, gsi::Covariance<MODEL>>
         makerGSI_Covariance_("gsi covariance");
  static SaberBlockMaker<MODEL, gsi::Interpolation<MODEL>>
         makerGSI_Interpolation_("gsi interpolation to model grid");

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
