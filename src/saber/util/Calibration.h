/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/Variables.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace util {

class calibrationWriteParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(calibrationWriteParameters, oops::Parameters)

 public:
  oops::RequiredParameter<std::string> covName{"covariance name",
    "covariance name global attribute", this};
  oops::RequiredParameter<std::string> mpiPattern{"mpi pattern",
    "mpi pattern", this};
  oops::RequiredParameter<std::string> filePath{"file path",
    "Path to written file", this};
};

void createCalibrationNetCDFHeaderInput(const eckit::LocalConfiguration & conf,
                                        const std::string & statsType,
                                        const std::string & binType,
                                        const oops::Variables & vars,  // variables with levels;
                                        const bool & doingCalibration,
                                        eckit::LocalConfiguration & netCDFConf);


}  // namespace util

