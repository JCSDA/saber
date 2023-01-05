/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/bump/type_bump_parameters.h"

#include "eckit/config/LocalConfiguration.h"

#include "saber/bump/BUMP.h"

namespace saber {
namespace bump {

// -----------------------------------------------------------------------------

void bump_config_init_f90(eckit::LocalConfiguration * config) {
  if (!config->empty()) {
    std::cout << "bump_config_init_f90: LocalConfiguration should be empty" << std::endl;
    std::abort();
  }

  BUMPParameters().serialize(*config);
}

// -----------------------------------------------------------------------------

}  // namespace bump
}  // namespace saber
