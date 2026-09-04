/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "saber/bifourier/BifourierTransformBase.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

class BifourierTransformStore {
 public:
  static const std::string classname()
    {return "saber::bifourier::BifourierTransformStore";}

  // Constructor
  BifourierTransformStore();

  // Destructor
  ~BifourierTransformStore();

  // Accessors

  // Store
  static std::vector<std::shared_ptr<BifourierTransformBase>> & transforms();

  // Return or create spectral transform from a grid-point function space (StructuredColumns)
  std::shared_ptr<BifourierTransformBase> setupTransform(
    const oops::GeometryData &,
    const oops::Variables &,
    const BifourierTransformParameters &) const;

  // Retrieve an existing spectral transform from a spectral function space (PointCloud)
  std::shared_ptr<BifourierTransformBase> retrieveTransform(
    const oops::GeometryData &) const;

  // Retrieve an existing spectral transform from a spectral function space (PointCloud),
  // for a specific number of variables/levels
  std::shared_ptr<BifourierTransformBase> retrieveTransform(
    const oops::GeometryData &,
    const oops::Variables &) const;
};

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
