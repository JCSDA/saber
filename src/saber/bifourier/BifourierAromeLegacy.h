/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#pragma once

#include <vector>

#include "atlas/field.h"

#include "eckit/mpi/Comm.h"

#include "saber/bifourier/BifourierTransform.h"

namespace saber {
namespace bifourier {

class BifourierVorToPbReadParameters;
class BifourierBalanceReadParameters;
class BifourierCovarianceReadParameters;

namespace arome_legacy {

// -----------------------------------------------------------------------------

void readVorToPb(const eckit::mpi::Comm &,
                 const BifourierVorToPbReadParameters &,
                 const BifourierTransform &,
                 std::vector<double> &);

// -----------------------------------------------------------------------------

void readBalance(const eckit::mpi::Comm &,
                 const oops::Variables &,
                 const BifourierBalanceReadParameters &,
                 const BifourierTransform &,
                 atlas::FieldSet &);

// -----------------------------------------------------------------------------

void readCovariance(const eckit::mpi::Comm &,
                    const oops::Variables &,
                    const BifourierCovarianceReadParameters &,
                    const BifourierTransform &,
                    atlas::FieldSet &);

// -----------------------------------------------------------------------------

}  // namespace arome_legacy
}  // namespace bifourier
}  // namespace saber
