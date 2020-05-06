/*
 * (C) Copyright 2020- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_INTERPOLATION_INTERPOLATORFACTORY_H_
#define SABER_INTERPOLATION_INTERPOLATORFACTORY_H_

#include<memory>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "eckit/config/LocalConfiguration.h"
#include "oops/generic/Interpolator.h"
#include "oops/generic/InterpolatorFactoryBase.h"

namespace saber {

// -----------------------------------------------------------------------------

/*! \brief Factory for creating oops::Interpolator objects
 *
 * This extends the oops::InterpolatorFactory to include bump interpolation
 * as an additional option
 *
 *\date May, 2020: initial implementation, M. Miesch (JCSDA)
 *
 */

class InterpolatorFactory : public oops::InterpolatorFactoryBase {
 public:
  std::string classname() const override {return "saber::InterpolatorFactory";}

  InterpolatorFactory() { }
  ~InterpolatorFactory() { }

  std::unique_ptr<oops::Interpolator> create(eckit::LocalConfiguration &,
                              const atlas::FunctionSpace &,
                              const atlas::FunctionSpace &,
                              const atlas::field::FieldSetImpl * = nullptr)
                              const override;
};

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_INTERPOLATION_INTERPOLATORFACTORY_H_
