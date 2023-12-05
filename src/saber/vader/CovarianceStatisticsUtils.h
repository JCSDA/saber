/*
 * (C) Crown Copyright 2022 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#pragma once

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "oops/base/Variables.h"

#include "saber/vader/PressureParameters.h"

namespace saber {
namespace vader {

atlas::Field createGpRegressionMatrices(const std::string &,
                                        const std::size_t,
                                        const std::size_t);

std::vector<double> interpWeights(std::vector<std::vector<double>> &,
                                  std::vector<std::vector<double>> &,
                                  double);

atlas::Field createGpRegressionWeights(const atlas::FunctionSpace &,
                                       const atlas::FieldSet &,
                                       const std::string &,
                                       const std::size_t,
                                       const std::size_t);

void interpMuStats(atlas::FieldSet &,
                   const atlas::Field &);

void populateMuA(atlas::FieldSet &,
                 const atlas::Field &);

atlas::FieldSet createGpRegressionStats(const atlas::FunctionSpace &,
                                        const atlas::FieldSet &,
                                        const oops::Variables &,
                                        const GpToHpCovarianceParameters &);
}  // namespace vader
}  // namespace saber
