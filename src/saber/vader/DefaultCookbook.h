/*
 * (C) Copyright 2024- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "vader/recipes/AirPressureExtendedUpByOne.h"
#include "vader/recipes/AirPressureThickness.h"
#include "vader/recipes/AirTemperature.h"
#include "vader/recipes/AirVirtualTemperature.h"
#include "vader/recipes/CloudIceMixingRatio.h"
#include "vader/recipes/CloudLiquidMixingRatio.h"
#include "vader/recipes/DryAirDensity.h"
#include "vader/recipes/DryAirDensityLevelsMinusOne.h"
#include "vader/recipes/HydrostaticExnerLevels.h"
#include "vader/recipes/HydrostaticPressureLevels.h"
#include "vader/recipes/LogDerivativeSaturationVaporPressure.h"
#include "vader/recipes/RainMixingRatio.h"
#include "vader/recipes/SaturationSpecificHumidity.h"
#include "vader/recipes/SaturationVaporPressure.h"
#include "vader/recipes/TotalMixingRatio.h"
#include "vader/recipes/TotalRelativeHumidity.h"
#include "vader/recipes/TotalWater.h"
#include "vader/recipes/VirtualPotentialTemperature.h"
#include "vader/recipes/WaterVaporMixingRatioWrtMoistAirAndCondensedWater.h"
#include "vader/vader.h"

namespace saber {

static const vader::cookbookConfigType saberDefaultCookbook = {
        {oops::Variable{"air_temperature"},
                                   {vader::AirTemperature_A::Name}},
        {oops::Variable{"air_pressure_levels"},
                                   {vader::AirPressureExtendedUpByOne_A::Name}},
        {oops::Variable{"mass_content_of_cloud_ice_in_atmosphere_layer"},
                                   {vader::CloudIceMixingRatio_A::Name}},
        {oops::Variable{"mass_content_of_cloud_liquid_water_in_atmosphere_layer"},
                                   {vader::CloudLiquidMixingRatio_A::Name}},
        {oops::Variable{"dry_air_density"},
                                   {vader::DryAirDensity_A::Name}},
        {oops::Variable{"dry_air_density_levels_minus_one"},
                                   {vader::DryAirDensityLevelsMinusOne_A::Name}},
        {oops::Variable{"dlsvpdT"},
                                   {vader::LogDerivativeSaturationVaporPressure_A::Name}},
        {oops::Variable{"hydrostatic_exner_levels"},
                                   {vader::HydrostaticExnerLevels_A::Name}},
        {oops::Variable{"hydrostatic_pressure_levels"},
                                   {vader::HydrostaticPressureLevels_A::Name}},
        {oops::Variable{"qrain"},  {vader::RainMixingRatio_A::Name}},
        {oops::Variable{"qsat"},   {vader::SaturationSpecificHumidity_A::Name}},
        {oops::Variable{"m_t"},    {vader::TotalMixingRatio_A::Name}},
        {oops::Variable{"qt"},     {vader::TotalWater_A::Name}},
        {oops::Variable{"rht"},    {vader::TotalRelativeHumidity_A::Name}},
        {oops::Variable{"specific_humidity"},
                       {vader::WaterVaporMixingRatioWrtMoistAirAndCondensedWater_A::Name}},
        {oops::Variable{"svp"},    {vader::SaturationVaporPressure_A::Name}},
        {oops::Variable{"virtual_potential_temperature"},
                                   {vader::VirtualPotentialTemperature_B::Name}},
        {oops::Variable{"virtual_temperature"},
                                   {vader::AirVirtualTemperature_A::Name}},
    };

}  // namespace saber
