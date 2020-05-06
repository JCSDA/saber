/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_INTERPOLATION_INSTANTIATEINTERPOLATORFACTORY_H_
#define SABER_INTERPOLATION_INSTANTIATEINTERPOLATORFACTORY_H_

#include "oops/generic/instantiateInterpolatorFactory.h"
#include "oops/base/InterpolatorBase.h"
#include "saber/oops/InterpolatorBump.h"

namespace saber {

void instantiateInterpolatorFactory() {
  oops::instantiateInterpolatorFactory();
  static oops::InterpolatorMaker<InterpolatorBump>
                        makerBumpInterpolator_("bump");
}

}  // namespace saber

#endif  // SABER_INTERPOLATION_INSTANTIATEINTERPOLATORFACTORY_H_
