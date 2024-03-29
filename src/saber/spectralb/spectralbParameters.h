/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_SPECTRALB_SPECTRALBPARAMETERS_H_
#define SABER_SPECTRALB_SPECTRALBPARAMETERS_H_

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
template <typename MODEL> class spectralbParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(spectralbParameters, oops::Parameters)

 public:
  // required
  oops::RequiredParameter<std::string> covarianceFile{"covariance_file",
  "Covariance file", this};
  oops::RequiredParameter<std::string> gaussGridUid{"gauss_grid_uid",
  "Gauss Grid UID", this};
  oops::RequiredParameter<std::vector<std::string>> umatrixNetCDFNames{
    "umatrix_netcdf_names", "umatrix netcdf names", this};

  // optional
  // TO DO - refactor SABER so that the two options below could be
  //         done both efficiently and flexibly within separate SABER blocks.
  //         This will probably require SABER blocks to have common initialisation
  //         blocks so that initialisation/setup steps are not duplicated.
  oops::OptionalParameter<std::string> outputGridUid{"output_grid_uid",
  "Output Grid UID", this};
  oops::OptionalParameter<bool> varianceOpt{"variance_opt",
  "variance option", this};
};
// -----------------------------------------------------------------------------
}  // namespace saber

#endif  // SABER_SPECTRALB_SPECTRALBPARAMETERS_H_
