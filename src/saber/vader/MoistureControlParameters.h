/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_VADER_MOISTURECONTROLPARAMETERS_H_
#define SABER_VADER_MOISTURECONTROLPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/base/ModelSpaceCovarianceParametersBase.h"
#include "oops/base/ParameterTraitsVariables.h"

#include "oops/util/DateTime.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace saber {

// -----------------------------------------------------------------------------
/// \brief Parameters passed to the Error Covariance class.
template <typename MODEL> class moisturecontrolParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(moisturecontrolParameters, oops::Parameters)

 public:
  // required
  oops::RequiredParameter<std::string> covariance_file_path{"covariance file path", this};
  oops::Parameter<int> mu_bins{"rht bins", "relative humidity bins", 30, this};
};
// -----------------------------------------------------------------------------
}  // namespace saber

#endif  // SABER_VADER_MOISTURECONTROLPARAMETERS_H_
