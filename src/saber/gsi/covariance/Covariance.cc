/*
 * (C) Copyright 2022 United States Government as represented by the Administrator of the National
 *     Aeronautics and Space Administration
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/gsi/covariance/Covariance.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/library.h"
#include "atlas/runtime/Log.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/Variables.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "saber/gsi/covariance/Covariance.interface.h"
#include "saber/gsi/grid/Grid.h"

namespace saber {
namespace gsi {

// -------------------------------------------------------------------------------------------------

static SaberCentralBlockMaker<Covariance> makerCovariance_("gsi covariance");

// -------------------------------------------------------------------------------------------------

Covariance::Covariance(const oops::GeometryData & geometryData,
                       const oops::Variables & centralVars,
                       const eckit::Configuration & covarConf,
                       const Parameters_ & params,
                       const oops::FieldSet3D & xb,
                       const oops::FieldSet3D & fg)
  : SaberCentralBlockBase(params, xb.validTime())
{
}

// -------------------------------------------------------------------------------------------------

Covariance::~Covariance() {
}

// -------------------------------------------------------------------------------------------------

}  // namespace gsi
}  // namespace saber
