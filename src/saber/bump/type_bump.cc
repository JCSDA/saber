/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/bump/type_bump.h"

#include "eckit/config/Configuration.h"

#include "saber/bump/BUMP.h"

namespace saber {
namespace bump {

void bump_config_init_f90(eckit::LocalConfiguration * config) {
  BUMPParameters().serialize(*config);
}

}  // namespace bump
}  // namespace saber
