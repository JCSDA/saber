/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/BifourierTransformStore.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

#include "saber/bifourier/BifourierUtilities.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

static std::vector<std::shared_ptr<BifourierTransform>> transformsVector;

// -----------------------------------------------------------------------------

std::vector<std::shared_ptr<BifourierTransform>> & BifourierTransformStore::transforms() {
  return transformsVector;
}

// -----------------------------------------------------------------------------

std::shared_ptr<BifourierTransform> BifourierTransformStore::setupTransform(
  const oops::GeometryData & gdata,
  const oops::Variables & vars,
  const eckit::Configuration & conf) const {
  oops::Log::trace() << classname() << "::setupTransform starting" << std::endl;

  // Check function space type
  ASSERT(gdata.functionSpace().type() == "StructuredColumns");

  // Get grid UID
  const atlas::functionspace::StructuredColumns fs(gdata.functionSpace());
  const std::string gridUid = fs.grid().uid();

  // Compare with existing grid UIDs
  for (auto it = transforms().begin(); it != transforms().end(); ++it) {
    if ((*it)->gridUid() == gridUid) {
      oops::Log::info() << "Info     : Retrieved Bifourier transform with grid UID: " << gridUid
        << std::endl;
      oops::Log::trace() << classname() << "::setupTransform done" << std::endl;
      return (*it);
    }
  }

  // Create transform
  std::shared_ptr<BifourierTransform> transform;
  transform.reset(new BifourierTransform(gdata, gridUid, vars, conf));

  // Insert new transform
  transforms().push_back(transform);

  oops::Log::trace() << classname() << "::setupTransform done" << std::endl;
  return *(std::prev(transforms().end()));
}

// -----------------------------------------------------------------------------

std::shared_ptr<BifourierTransform> BifourierTransformStore::retrieveTransform(
  const oops::GeometryData & gdata) const {
  oops::Log::trace() << classname() << "::retrieveTransform starting" << std::endl;

  // Check function space type
  ASSERT(gdata.functionSpace().type() == "PointCloud");

  // Generate spectral UID
  const std::string specUid = generateSpectralUid(gdata.functionSpace(), gdata.comm());

  // Compare with existing spectral UIDs
  for (auto it = transforms().begin(); it != transforms().end(); ++it) {
    if ((*it)->specUid() == specUid) {
      oops::Log::info() << "Info     : Retrieved Bifourier transform with spectral UID: " << specUid
        << std::endl;

      oops::Log::trace() << classname() << "::setupTransform done" << std::endl;
      return (*it);
    }
  }

  // Cannot find existing spectral space
  throw eckit::Exception("cannot find existing spectral space", Here());

  oops::Log::trace() << classname() << "::retrieveTransform done" << std::endl;
  return *(std::prev(transforms().end()));
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
