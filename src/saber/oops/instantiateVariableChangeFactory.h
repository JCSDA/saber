/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef SABER_OOPS_INSTANTIATEVARIABLECHANGEFACTORY_H_
#define SABER_OOPS_INSTANTIATEVARIABLECHANGEFACTORY_H_

#include "oops/generic/instantiateVariableChangeFactory.h"

#include "saber/oops/PsiChiToUVVariableChange.h"
#include "saber/oops/StatsVariableChange.h"
#include "saber/oops/StdDevVariableChange.h"

namespace saber {

template <typename MODEL>
void instantiateVariableChangeFactory() {
  oops::instantiateVariableChangeFactory<MODEL>();
  static oops::LinearVariableChangeMaker<MODEL, StatsVariableChange<MODEL> >
                        makerStatsVarChange_("StatsVariableChange");
  static oops::LinearVariableChangeMaker<MODEL, StdDevVariableChange<MODEL> >
                        makerStdDev_("StdDev");
  static oops::LinearVariableChangeMaker<MODEL, PsiChiToUVVariableChange<MODEL> >
                        makerPsiChiToUV_("PsiChiToUV");
}

}  // namespace saber

#endif  // SABER_OOPS_INSTANTIATEVARIABLECHANGEFACTORY_H_
