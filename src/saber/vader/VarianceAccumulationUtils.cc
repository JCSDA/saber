/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/VarianceAccumulationUtils.h"

#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/meshgenerator.h"
#include "atlas/option.h"
#include "atlas/util/Earth.h"
#include "atlas/util/GaussianLatitudes.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/base/FieldSet3D.h"
#include "oops/util/AtlasArrayUtil.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"


namespace saber {
namespace vader {

using atlas::array::make_view;
using atlas::idx_t;

atlas::Field CSCellArea(const eckit::mpi::Comm & localComm,
                        const atlas::functionspace::NodeColumns & nodeFspace,
                        const std::size_t & sizeOwned,
                        const std::string& fieldName)
{
  constexpr auto degreesToRadians = M_PI / 180.;
  const auto cubedSphereGrid = atlas::CubedSphereGrid(nodeFspace.mesh().grid());
  const auto lonlat = make_view<double, 2>(nodeFspace.lonlat());
  constexpr double tolerance = 1e-8;

  const auto& proj = cubedSphereGrid.cubedSphereProjection();

  atlas::Field gridAreaField =
    atlas::Field(fieldName,  atlas::array::make_datatype<double>(),
                             atlas::array::make_shape(1, sizeOwned));
  auto gridAreaView = make_view<double, 2>(gridAreaField);

  // Area is calculated by transforming the grid area of a cell on a cube onto the
  // Mercator projection and then converting it onto the sphere by applying the cosine.
  // This is only an approximation to the cell area and only converges linearly as a function
  // of resolution.
  double totalArea(0.0);
  for (std::size_t i = 0; i < gridAreaView.size(); i++) {
    const atlas::PointLonLat ll = atlas::PointLonLat(lonlat(i, atlas::LON),
                                                     lonlat(i, atlas::LAT));
    if (std::abs(std::abs(ll.lat()) - 90.0) < tolerance) {
      gridAreaView(0, i) = 1.0;
    } else {
      const double cosLat = std::cos(degreesToRadians * ll.lat());
      // The absolute value of the Jacobian determinant at p gives us
      // the factor by which the function f expands or shrinks volumes near a
      // given point. We are looking at the reciprocal here as we are
      // considering the inverse of the projection.
      const double gridJacDet = 1.0/proj.jacobian(ll).determinant();
      gridAreaView(0, i) = cosLat * gridJacDet;
    }
    totalArea += gridAreaView(0, i);
  }
  localComm.allReduceInPlace(totalArea, eckit::mpi::sum());
  for (std::size_t i = 0; i < sizeOwned; i++) {
    gridAreaView(0, i) =   gridAreaView(0, i) / totalArea;
  }
  return gridAreaField;
}

atlas::Field GaussCellArea(const eckit::mpi::Comm & localComm,
                           const atlas::functionspace::StructuredColumns & structCols,
                           const std::size_t & sizeOwned,
                           const std::string& fieldName)
{
  const auto structGrid = atlas::StructuredGrid(structCols.grid());
  auto lonlatView = make_view<double, 2>(structCols.lonlat());
  auto northSouthPts = structGrid.ny();
  std::vector<double> latitudes(northSouthPts, 0.0);
  std::vector<double> weights(northSouthPts, 0.0);
  // The function below calculates the Gaussian latitudes (latitudes) and a
  // normalized surface area assigned to each Gaussian latitude (weights)
  // The normalized surface area weights sum to 1.
  atlas::util::gaussian_quadrature_npole_spole(northSouthPts/2,
                                               latitudes.data(),
                                               weights.data());
  atlas::Field gridAreaField =
    atlas::Field(fieldName,  atlas::array::make_datatype<double>(),
                             atlas::array::make_shape(1, sizeOwned));
  auto gridAreaView = make_view<double, 2>(gridAreaField);
  gridAreaView.assign(0.0);

  double totalArea(0.0);
  for (std::size_t i = 0; i < sizeOwned; i++) {
    auto ll = atlas::PointLonLat(lonlatView(i, atlas::LON), lonlatView(i, atlas::LAT));

    // find index of lat
    for (std::size_t j = 0;  j < latitudes.size(); j++) {
      if (std::abs(latitudes[j] - ll.lat()) < __FLT_EPSILON__) {
        gridAreaView(0, i) = weights[j]  / static_cast<double>(structGrid.nx(j));
        break;
      }
    }
    totalArea += gridAreaView(0, i);
  }
  localComm.allReduceInPlace(totalArea, eckit::mpi::sum());
  for (std::size_t i = 0; i < sizeOwned; i++) {
    gridAreaView(0, i) =   gridAreaView(0, i) / totalArea;
  }
  return gridAreaField;
}

void computeBinnedVariances(const eckit::mpi::Comm & localComm,
                            const std::string & tag,
                            const atlas::Field & gridptAreaWeights,
                            const idx_t & sizeOwned,
                            const atlas::FieldSet & fset,
                            const std::size_t & localRoot,
                            const std::size_t & nbins,
                            atlas::FieldSet & variances)
{
  ASSERT(nbins == 1);
  ASSERT(localRoot == 0);
  auto gridAreaView = make_view<double, 2>(gridptAreaWeights);

  for (const auto & field : fset) {
    auto fldView = make_view<double, 2>(fset[field.name()]);
    auto fldVariance = atlas::Field(field.name(), atlas::array::make_datatype<double>(),
      atlas::array::make_shape(nbins, fldView.shape(1)));
    auto fldVarianceView = make_view<double, 2>(fldVariance);
    fldVarianceView.assign(0.0);

    for (idx_t jl = 0; jl < fldView.shape(1); ++jl) {
      for (idx_t jn = 0; jn < sizeOwned; ++jn) {
        fldVarianceView(0, jl) += fldView(jn, jl) * fldView(jn, jl) * gridAreaView(0, jn);
      }
    }

    // Note that the variances will only have the correct value on the root PE.
    util::gatherSum(localComm, localRoot, fldVarianceView);

    variances.add(fldVariance);
  }
}

// This should work for both cubed-sphere and gauss meshes
// TODO(Marek) - replace this code with calibration code
void printInstantBinnedVariances(const eckit::mpi::Comm & globalComm,
                                 const std::size_t & globalRoot,
                                 const std::string & tag,
                                 const std::size_t & nbins,
                                 const atlas::FieldSet & variances)
{
  ASSERT(nbins == 1);
  ASSERT(globalRoot == 0);
  // print variance for each field
  for (const auto & f : variances) {
    auto fldVarianceView = make_view<double, 2>(f);
    std::vector<double> variance(nbins, 0.0);
    for (std::size_t b = 0; b < nbins; ++b) {
      for (idx_t l = 0; l < fldVarianceView.shape(1); ++l) {
        variance[b] += fldVarianceView(b, l);
      }
      variance[b] /= static_cast<double>(fldVarianceView.shape(1));
    }

    for (std::size_t b = 0; b < nbins; ++b) {
      oops::Log::test() << tag
                        << " ; name = " << f.name()
                        << " ; bin index = " << b
                        << " ; binned variance = " << variance[b]
                        << std::endl;
    }
  }
}

// TODO(Marek) - replace this code with calibration code
void computeVarianceFieldSetInstant(const eckit::mpi::Comm & comm,
                                    const std::string & tag,
                                    const atlas::FieldSet & fset,
                                    const std::size_t & root,
                                    const std::size_t & nbins,
                                    atlas::FieldSet & variances) {
  ASSERT(nbins == 1);
  ASSERT(root == 0);
  bool computeVariance(false);
  if (fset[0].functionspace().type() == "StructuredColumns") {
    auto scFspace = atlas::functionspace::StructuredColumns(fset[0].functionspace());
    std::string s = scFspace.grid().name().substr(0, 1);
    // "F" denotes a regular Gaussian grid.
    if (s == "F") {
      atlas::Field gridptAreaWeights = GaussCellArea(comm,
                                                     scFspace,
                                                     scFspace.sizeOwned(),
                                                     "grid point area weights");
      computeBinnedVariances(comm,
                             tag,
                             gridptAreaWeights,
                             static_cast<idx_t>(scFspace.sizeOwned()),
                             fset,
                             root,
                             nbins,
                             variances);
      computeVariance = true;
    }
  } else if (fset[0].functionspace().type() == "NodeColumns") {
    std::string s = atlas::functionspace::NodeColumns(
      fset[0].functionspace()).mesh().grid().name();
    if (s.substr(0, 6) == "CS-LFR") {
      auto csFspace =
        atlas::functionspace::CubedSphereNodeColumns(fset[0].functionspace());
      atlas::Field gridptAreaWeights =
        CSCellArea(comm, fset[0].functionspace(),
                   static_cast<const std::size_t>(csFspace.sizeOwned()),
                   std::string{"grid point area weights"});
      computeBinnedVariances(comm,
                             tag,
                             gridptAreaWeights,
                             static_cast<idx_t>(csFspace.sizeOwned()),
                             fset,
                             root,
                             nbins,
                             variances);
      computeVariance = true;
    }
  }
  if (!computeVariance) {
    oops::Log::error() << "WriteVariances:: variances calculations not supported "
                       << "for this specific mesh" << std::endl;
  }
}

std::size_t updateVerticalCovariances(const eckit::LocalConfiguration & netCDFConf,
                                      const atlas::FieldSet & binningData,
                                      const oops::FieldSets & ensFieldSet,
                                      const std::size_t priorSampleSize,
                                      atlas::FieldSet & ensembleStats) {
  if (priorSampleSize > 0) {
    // assuming read done before and multiplying by number of prior samples
    for (atlas::Field & vertcov : ensembleStats) {
      if (vertcov.name().compare(0, 19, "vertical covariance") == 0) {
        auto vertCovView = make_view<double, 3>(vertcov);
        for (idx_t i = 0; i < vertcov.shape()[0]; ++i) {
          for (idx_t j = 0; j < vertcov.shape()[1]; ++j) {
            for (idx_t k = 0; k < vertcov.shape()[2]; ++k) {
              vertCovView(i, j, k) *= static_cast<double>(priorSampleSize);
            }
          }
        }
      }
    }
  }

  for (atlas::Field & vertcov : ensembleStats) {
    if (vertcov.name().compare(0, 19, "vertical covariance") == 0) {
      const std::string var1 =
        util::getAttributeValue<std::string>(netCDFConf, vertcov.name(),
                                             "variable name 1");
      const std::string var2 =
        util::getAttributeValue<std::string>(netCDFConf, vertcov.name(),
                                             "variable name 2");
      std::string bintype =
        util::getAttributeValue<std::string>(netCDFConf, vertcov.name(),
                                             "binning type");

      std::string binindx = bintype + " local indices";
      auto wgtView = make_view<double, 2>(binningData[bintype + " weights"]);
      auto indxView = make_view<std::int32_t, 2>(binningData[binindx]);
      auto vertCovView = make_view<double, 3>(vertcov);

      for (size_t jj = 0; jj < ensFieldSet.size(); ++jj) {
        const auto & fs = ensFieldSet[jj];
        auto var1View = make_view<const double, 2>(fs[var1]);
        auto var2View = make_view<const double, 2>(fs[var2]);

        for (idx_t b = 0; b < vertcov.shape()[0]; ++b) {
          for (idx_t i = 0; i < binningData[binindx].shape()[1]; ++i) {
            for (idx_t j = 0; j < vertcov.shape()[1]; ++j) {
              for (idx_t k = 0; k < vertcov.shape()[2]; ++k) {
                std::int32_t l = indxView(b, i);
                vertCovView(b, j, k) += var1View(l, j) * var2View(l, k)
                                      * wgtView(b, i);
              }
            }
          }
        }
      }
    }
  }

  const std::size_t updatedSampleSize = priorSampleSize + ensFieldSet.size();

  const double recipPriorSampleSize = 1.0/static_cast<double>(updatedSampleSize);
  for (atlas::Field & vertcov : ensembleStats) {
    if (vertcov.name().compare(0, 19, "vertical covariance") == 0) {
      auto vertCovView = make_view<double, 3>(vertcov);
      for (idx_t i = 0; i < vertcov.shape()[0]; ++i) {
        for (idx_t j = 0; j < vertcov.shape()[1]; ++j) {
          for (idx_t k = 0; k < vertcov.shape()[2]; ++k) {
            vertCovView(i, j, k) *= recipPriorSampleSize;
          }
        }
      }
    }
  }
  return updatedSampleSize;
}

std::size_t updateVariances(const eckit::LocalConfiguration & netCDFConf,
                     const atlas::FieldSet & binningData,
                     const oops::FieldSets & ensFieldSet,
                     const std::size_t priorSampleSize,
                     atlas::FieldSet & ensembleStats) {
  if (priorSampleSize > 0) {
    // assuming read done before and multiplying by number of prior samples
    for (atlas::Field & variance : ensembleStats) {
      if (variance.name().compare(0, 8, "variance") == 0) {
        auto varianceView = make_view<double, 2>(variance);
        for (idx_t i = 0; i < variance.shape()[0]; ++i) {
          for (idx_t j = 0; j < variance.shape()[1]; ++j) {
            varianceView(i, j) *= static_cast<double>(priorSampleSize);
          }
        }
      }
    }
  }

  for (atlas::Field & variance : ensembleStats) {
    if (variance.name().compare(0, 8, "variance") == 0) {
      const std::string var1 =
        util::getAttributeValue<std::string>(netCDFConf, variance.name(),
                                             "variable name 1");
      const std::string var2 =
        util::getAttributeValue<std::string>(netCDFConf, variance.name(),
                                             "variable name 2");
      std::string bintype =
        util::getAttributeValue<std::string>(netCDFConf, variance.name(),
                                             "binning type");
      std::string binindx = bintype + " local indices";
      auto wgtView = make_view<double, 2>(binningData[bintype + " weights"]);
      auto indxView = make_view<std::int32_t, 2>(binningData[binindx]);
      auto varianceView = make_view<double, 2>(variance);

      for (size_t jj = 0; jj < ensFieldSet.size(); ++jj) {
        const auto & fs = ensFieldSet[jj];
        auto var1View = make_view<const double, 2>(fs[var1]);
        auto var2View = make_view<const double, 2>(fs[var2]);

        for (idx_t b = 0; b < variance.shape()[0]; ++b) {
          for (idx_t i = 0; i < binningData[binindx].shape()[1]; ++i) {
            for (idx_t j = 0; j < variance.shape()[1]; ++j) {
              std::int32_t l = indxView(b, i);
              varianceView(b, j) += var1View(l, j) * var2View(l, j)
                                    * wgtView(b, i);
            }
          }
        }
      }
    }
  }

  const std::size_t updatedSampleSize = priorSampleSize + ensFieldSet.size();

  const double recipPriorSampleSize = 1.0/static_cast<double>(updatedSampleSize);

  for (atlas::Field & variance : ensembleStats) {
    if (variance.name().compare(0, 8, "variance") == 0) {
      auto varianceView = make_view<double, 2>(variance);
      for (idx_t i = 0; i < variance.shape()[0]; ++i) {
        for (idx_t j = 0; j < variance.shape()[1]; ++j) {
          varianceView(i, j) *= recipPriorSampleSize;
        }
      }
    }
  }
  return updatedSampleSize;
}

void copyFieldSet(const atlas::FieldSet & otherFset,
                  atlas::FieldSet & fset) {
  ASSERT(fset.empty());
  for (idx_t ivar = 0; ivar < otherFset.size(); ++ivar) {
    if (otherFset[ivar].rank() == 2) {
      auto fld =
        atlas::Field(otherFset.field_names()[ivar],
                     atlas::array::make_datatype<double>(),
                     atlas::array::make_shape(otherFset[ivar].shape(0),
                                              otherFset[ivar].shape(1)));
      auto fView = make_view<double, 2>(fld);
      fView.assign(0.0);
      auto otherView = make_view<const double, 2>(otherFset[fld.name()]);
      for (idx_t jn = 0; jn < fld.shape(0); ++jn) {
        for (idx_t jl = 0; jl < fld.shape(1); ++jl) {
          fView(jn, jl) = otherView(jn, jl);
        }
      }
      fld.metadata() = otherFset[fld.name()].metadata();
      fset.add(fld);
    } else if (otherFset[ivar].rank() == 3) {
      auto fld =
        atlas::Field(otherFset.field_names()[ivar],
                     atlas::array::make_datatype<double>(),
                     atlas::array::make_shape(otherFset[ivar].shape(0),
                                              otherFset[ivar].shape(1),
                                              otherFset[ivar].shape(2)));
      auto fView = make_view<double, 3>(fld);
      fView.assign(0.0);
      auto otherView = make_view<const double, 3>(otherFset[fld.name()]);
      for (idx_t jn = 0; jn < fld.shape(0); ++jn) {
        for (idx_t jl = 0; jl < fld.shape(1); ++jl) {
          for (idx_t jl2 = 0; jl2 < fld.shape(2); ++jl2) {
            fView(jn, jl, jl2) = otherView(jn, jl, jl2);
          }
        }
      }
      fld.metadata() = otherFset[fld.name()].metadata();
      fset.add(fld);
    } else if (otherFset[ivar].rank() == 1) {
      auto fld =
        atlas::Field(otherFset.field_names()[ivar],
                     atlas::array::make_datatype<double>(),
                     atlas::array::make_shape(otherFset[ivar].shape(0)));
      auto fView = make_view<double, 1>(fld);
      fView.assign(0.0);
      auto otherView = make_view<const double, 1>(otherFset[fld.name()]);
      for (idx_t jn = 0; jn < fld.shape(0); ++jn) {
        fView(jn) = otherView(jn);
      }
      fld.metadata() = otherFset[fld.name()].metadata();
      fset.add(fld);
    } else {
      throw eckit::UserError("copy FieldSet works only with 1D and 2D fields", Here());
    }
  }
}

}  // namespace vader
}  // namespace saber

