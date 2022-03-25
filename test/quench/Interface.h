/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QUENCH_INTERFACE_H_
#define QUENCH_INTERFACE_H_

#include <memory>

#include "atlas/functionspace.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace quench {

extern "C" {
  void fields_write_structuredcolumns_f90(const eckit::Configuration &, const oops::Variables &,
                                          atlas::field::FieldSetImpl *,
                                          atlas::field::FieldSetImpl *);
}

}  // namespace quench
#endif  // QUENCH_INTERFACE_H_
