/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "oops/base/FieldSets.h"

namespace saber {
namespace vader {

atlas::Field CSCellArea(const eckit::mpi::Comm & localComm,
                        const atlas::functionspace::NodeColumns & nodeFspace,
                        const std::size_t & sizeOwned,
                        const std::string& fieldName);

atlas::Field GaussCellArea(
    const eckit::mpi::Comm & localComm,
    const atlas::functionspace::StructuredColumns & structCols,
    const std::size_t & sizeOwned,
    const std::string& fieldName);

void printInstantBinnedVariances(const eckit::mpi::Comm & globalComm,
                                 const std::size_t & globalRoot,
                                 const std::string & tag,
                                 const std::size_t & nbins,
                                 const atlas::FieldSet & variances);

void computeVarianceFieldSetInstant(const eckit::mpi::Comm & comm,
                                    const std::string & tag,
                                    const atlas::FieldSet & fset,
                                    const std::size_t & root,
                                    const std::size_t & nbins,
                                    atlas::FieldSet & variances);

std::size_t updateVerticalCovariances(const eckit::LocalConfiguration & netCDFConf,
                                      const atlas::FieldSet & binningData,
                                      const oops::FieldSets & ensFieldSet,
                                      const std::size_t priorSampleSize,
                                      atlas::FieldSet & verticalCovariances);

std::size_t updateVariances(const eckit::LocalConfiguration & netCDFconf,
                            const atlas::FieldSet & binningData,
                            const oops::FieldSets & ensFieldSet,
                            const std::size_t priorSampleSize,
                            atlas::FieldSet & variances);

void copyFieldSet(const atlas::FieldSet & otherFset,
                  atlas::FieldSet & fset);
}  // namespace vader
}  // namespace saber
