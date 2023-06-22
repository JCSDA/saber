/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include <ostream>
#include <string>
#include <vector>

#include "atlas/util/Config.h"

#include "saber/quench/src/ModelData.h"

// -------------------------------------------------------------------------------------------------

namespace quench {

// -------------------------------------------------------------------------------------------------

ModelData::ModelData(const Geometry & geometry) {}

// -------------------------------------------------------------------------------------------------

ModelData::~ModelData() {}

// -------------------------------------------------------------------------------------------------

const eckit::LocalConfiguration ModelData::modelData() const {
  eckit::LocalConfiguration quenchModelData;
  quenchModelData.set("epsilon", 0.621957535);
  return quenchModelData;
}

// -------------------------------------------------------------------------------------------------

void ModelData::print(std::ostream & os) const {
  os << "quench::ModelData::modelData(): " << modelData();
}

// -------------------------------------------------------------------------------------------------

}  // namespace quench
