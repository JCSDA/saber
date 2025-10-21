/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "saber/bifourier/BifourierTransform.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

class BifourierTransformStore {
 public:
  static const std::string classname()
    {return "saber::bifourier::BifourierTransform";}

  // Constructor
  BifourierTransformStore();

  // Destructor
  ~BifourierTransformStore();

  // Accessors

  // Store
  static std::vector<std::shared_ptr<BifourierTransform>> & transforms();

  // Return or create spectral transform from a grid-point function space (StructuredColumns)
  std::shared_ptr<BifourierTransform> setupTransform(
    const oops::GeometryData &,
    const oops::Variables &,
    const eckit::Configuration &) const;

  // Retrieve an existing spectral transform from a spectral function space (PointCloud)
  std::shared_ptr<BifourierTransform> retrieveTransform(
    const oops::GeometryData &) const;
};

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
