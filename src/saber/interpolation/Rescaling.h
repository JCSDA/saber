/*
 * (C) Crown Copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "atlas/functionspace.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/Variables.h"

#include "saber/interpolation/AtlasInterpWrapper.h"


namespace saber {
namespace interpolation {

/// \brief Rescaling layer to mitigate variance loss in interpolation
/// \details This class computes and applies multiplicative normalization coefficients
///          after interpolation in a covariance or correlation matrix. It only
///          works with the atlas::interpolation::AtlasInterpolationWrapper interpolation.
///          This could have been defined as a SaberOuterBlock, but would have required
///          recomputing an atlas interpolation and redistribution. This is best used as a
///          member in an interpolation class to reuse what is already computed here.
class Rescaling{
 public:
  /// \brief Constructor from covariance profiles files generated by ErrorCovarianceToolbox
  Rescaling(const eckit::mpi::Comm & comm,
            const eckit::LocalConfiguration & conf,
            const oops::Variables & vars,
            const atlas::FunctionSpace & innerFspace,
            const atlas::FunctionSpace & outerFspace,
            const saber::interpolation::AtlasInterpWrapper & interp);

  /// \brief Constructor from rescaling file generated by prior run
  Rescaling(const eckit::mpi::Comm & comm,
            const eckit::LocalConfiguration & conf,
            const oops::Variables & vars,
            const atlas::FunctionSpace & outerFspace);

  /// \brief Constructor for an empty rescaling, does nothing
  Rescaling() {}

  void execute(oops::FieldSet3D & fieldSet) const;

 private:
  /// Multiplicative coefficients, diagonal of D in C -> D C D^t
  const atlas::FieldSet multiplicativeRescalingCoeffs_;

  std::string classname() const {return "Rescaling";}
};

}  // namespace interpolation
}  // namespace saber