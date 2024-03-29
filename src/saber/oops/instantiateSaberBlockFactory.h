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

#if ENABLE_GSIBEC
  #include "saber/gsi/covariance/GSI_Covariance.h"
  #include "saber/gsi/interpolation/GSI_Interpolation.h"
#endif

#include "saber/oops/ID.h"
#include "saber/oops/StdDev.h"
#if atlas_TRANS_FOUND
  #include "saber/spectralb/SPCTRL_Cov.h"
  #include "saber/spectralb/SPNOINTERP_Cov.h"
#endif
#if ENABLE_VADER
  #if ENABLE_VADER_MO
    #include "saber/vader/AirTemperatureSaberBlock.h"
    #include "saber/vader/DryAirDensitySaberBlock.h"
    #include "saber/vader/HydroBalSaberBlock.h"
    #include "saber/vader/HydrostaticExnerSaberBlock.h"
    #include "saber/vader/MoistIncrOpSaberBlock.h"
    #include "saber/vader/MoistureControlSaberBlock.h"
  #endif
#endif

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
#if ENABLE_GSIBEC
  static SaberBlockMaker<MODEL, gsi::Covariance<MODEL>>
         makerGSI_Covariance_("gsi covariance");
  static SaberBlockMaker<MODEL, gsi::Interpolation<MODEL>>
         makerGSI_Interpolation_("gsi interpolation to model grid");
#endif

  // Identity
  static SaberBlockMaker<MODEL, ID<MODEL> >
         makerID_("ID");

#if atlas_TRANS_FOUND
  // Spectral B
  static SaberBlockMaker<MODEL, spectralb::SPCTRL_COV<MODEL> >
         makerSPCTRL_COV_("SPCTRL_COV");
  static SaberBlockMaker<MODEL, spectralb::SPNOINTERP_COV<MODEL> >
         makerSPNOINTERP_COV_("SPNOINTERP_COV");
#endif

#if ENABLE_VADER
  // VADER-based blocks
  #if ENABLE_VADER_MO
    // MO-specific VADER-based blocks
    static SaberBlockMaker<MODEL, AirTemperatureSaberBlock<MODEL> >
           makerAirTemperatureSaberBlock_("mo_air_temperature");
    static SaberBlockMaker<MODEL, DryAirDensitySaberBlock<MODEL> >
           makerDryAirDensitySaberBlock_("mo_dry_air_density");
    static SaberBlockMaker<MODEL, HydroBalSaberBlock<MODEL> >
           makerHydroBalSaberBlock_("mo_hydro_bal");
    static SaberBlockMaker<MODEL, HydrostaticExnerSaberBlock<MODEL> >
           makerHydrostaticExnerSaberBlock_("mo_hydrostatic_exner");
    static SaberBlockMaker<MODEL, MoistureControlSaberBlock<MODEL> >
           makerMoistureControlBlock_("mo_moisture_control");
    static SaberBlockMaker<MODEL, MoistIncrOpSaberBlock<MODEL> >
           makerMoistIncrOpSaberBlock_("mo_moistincrop");
  #endif
#endif

  // StdDev
  static SaberBlockMaker<MODEL, StdDev<MODEL> >
         makerStdDev_("StdDev");
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_INSTANTIATESABERBLOCKFACTORY_H_
