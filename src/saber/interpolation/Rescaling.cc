/*
 * (C) Crown Copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <netcdf.h>

#include <cmath>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>  // for std::pair
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/field/for_each.h"
#include "atlas/functionspace.h"
#include "atlas/util/CoordinateEnums.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/linalg/SparseMatrix.h"
#include "eckit/mpi/Comm.h"

#include "oops/base/Variables.h"
#include "oops/util/AtlasArrayUtil.h"
#include "oops/util/FieldSetHelpers.h"

#include "saber/interpolation/AtlasInterpWrapper.h"
#include "saber/interpolation/Rescaling.h"

namespace {

// -----------------------------------------------------------------------------

auto readCovarianceProfiles(const std::string & filePath,
                            const oops::Variables & activeVars) {
  oops::Log::trace() << "readCovarianceProfiles starting" << std::endl;
  // Setup
  std::vector<std::string> dimNames;
  std::vector<atlas::idx_t> dimSizes;
  std::vector<std::vector<std::string>> dimNamesForEveryVar;
  oops::Variables vars;
  eckit::LocalConfiguration netcdfMetaData;
  std::vector<int> netcdfGeneralIDs;
  std::vector<int> netcdfDimIDs;
  std::vector<int> netcdfVarIDs;
  std::vector<std::vector<int>> netcdfDimVarIDs;

  util::atlasArrayInquire(filePath,
                          dimNames,
                          dimSizes,
                          vars,
                          dimNamesForEveryVar,
                          netcdfMetaData,
                          netcdfGeneralIDs,
                          netcdfDimIDs,
                          netcdfVarIDs,
                          netcdfDimVarIDs);

  std::unordered_map<std::string, std::vector<double>> distances;
  atlas::FieldSet covariances;
  for (const auto & var : activeVars) {
    // Check variables are in file
    for (const auto ext : {"distances", "levels", "covariances"}) {
      const auto extVar = var.name() + "_1_" + ext;
      if (!vars.has(extVar)) {
        std::stringstream errorMsg;
        errorMsg << "Could not find variable " << extVar
                 << " in " << filePath << std::endl;
        throw eckit::UserError(errorMsg.str(), Here());
      }
    }

    // Check all levels are present in levels variable in file
    const auto levVarName = var.name() + "_1_levels";
    auto pos = vars.find(levVarName);
    ASSERT(netcdfDimVarIDs[pos].size() == 1);
    const int nLevs = dimSizes[netcdfDimVarIDs[pos][0]];
    ASSERT(nLevs == var.getLevels());
    auto levArray = atlas::array::Array::create<int>(nLevs);
    auto levView = atlas::array::make_view<int, 1>(*levArray);
    util::atlasArrayReadData(netcdfGeneralIDs,
                             netcdfVarIDs[pos],
                             levView);
    int lev(0);
    for (int jlev = 0; jlev < nLevs; jlev++) {
      ASSERT(levView(jlev) == lev++);
    }

    // Read vector of distances
    const auto distVarName = var.name() + "_1_distances";
    pos = vars.find(distVarName);
    ASSERT(netcdfDimVarIDs[pos].size() == 1);
    const int nDist = dimSizes[netcdfDimVarIDs[pos][0]];
    auto distArray = atlas::array::Array::create<double>(nDist);
    auto distView = atlas::array::make_view<double, 1>(*distArray);
    util::atlasArrayReadData(netcdfGeneralIDs,
                             netcdfVarIDs[pos],
                             distView);
    // Initialize and fill distance vector
    std::vector<double> distance;
    distance.reserve(nDist);
    for (int jnode = 0; jnode < nDist; jnode++) {
      distance.push_back(distView(jnode));
    }
    distances[var.name()] = distance;

    // Read field of covariances
    const auto covVarName = var.name() + "_1_covariances";
    pos = vars.find(covVarName);
    ASSERT(netcdfDimVarIDs[pos].size() == 2);
    ASSERT(dimSizes[netcdfDimVarIDs[pos][0]] == nDist);
    ASSERT(dimSizes[netcdfDimVarIDs[pos][1]] == nLevs);
    auto covField = atlas::Field(var.name(),
                                 atlas::array::make_datatype<double>(),
                                 atlas::array::make_shape(nDist, nLevs));
    auto covView = atlas::array::make_view<double, 2>(covField);
    util::atlasArrayReadData(netcdfGeneralIDs,
                             netcdfVarIDs[pos],
                             covView);
    covariances.add(covField);
  }

  int retval;
  if ((retval = nc_close(netcdfGeneralIDs[0]))) {
    throw eckit::Exception(nc_strerror(retval), Here());
  }

  oops::Log::trace() << "readCovarianceProfiles about to exit..." << std::endl;
  return std::make_pair(distances, covariances);
}

// -----------------------------------------------------------------------------
// lookUp function to interpolate covariance value at a given distance.
//
// Input parameters:
// - xx: 1D input distance vector, sorted by increasing order.
// - yView: 2D array view so that yy := yView(:, lev) is the output vector of
//          covariances associated to xx.
//          xx and yy define a piecewise linear function f so that f(xx) == yy.
// - lev: vertical level to define yy from yView.
// - x: input distance, with xx[0] <= x <= xx[xx.size()-1].
//
// Returns:
//  The interpolated covariance value y == f(x).
double lookUpCovarianceAtDistance(const std::vector<double> & xx,
                                  const atlas::array::ArrayView<const double, 2> & yView,
                                  const size_t lev,
                                  const double x) {
  if (x == xx[0]) return yView(0, lev);

  // Find position of x in xx
  const auto ind = std::lower_bound(xx.cbegin(), xx.cend(), x);
  if (ind == xx.cend()) {
    std::stringstream errorMsg;
    errorMsg << "Received distance " << x << "m, but maximum input distance"
             << " in profile is " << xx[xx.size()-1];
    throw eckit::UserError(errorMsg.str(), Here());
  }
  const auto pos = std::distance(xx.cbegin(), ind);

  // Perform linear interpolation
  const double a = (xx[pos] - x) / (xx[pos] - xx[pos-1]);
  const double y = a * yView(pos-1, lev) + (1-a) * yView(pos, lev);
  return y;
}

// -----------------------------------------------------------------------------

std::tuple<atlas::FieldSet, atlas::FieldSet>
computeRescalingCoeffs(
        const oops::Variables & activeVars,
        const atlas::FunctionSpace & innerFspace,
        const saber::interpolation::AtlasInterpWrapper & interp,
        const std::unordered_map<std::string, std::vector<double>> & binnedDistances,
        const atlas::FieldSet & covariances,
        const std::vector<double> & alphas) {
  oops::Log::trace() << "computeRescalingCoeffs starting" << std::endl;

  const auto & interpMatrix = interp.getInterpolationMatrix();
  // matData is the array of all non-zero matrix coefficients
  // matInner is the array of column indexes of all non-zero matrix coefficients
  // matOuter is the array of indexes of matData corresponding to the start of each row
  const auto & matData = interpMatrix.data();
  const auto & matInner = interpMatrix.inner();
  const auto & matOuter = interpMatrix.outer();

  const auto & matchingOuterFspace = interp.getIntermediateFunctionSpace();

  atlas::FieldSet multFset;
  atlas::FieldSet addFset;
  int ivar = 0;
  for (const auto & var : activeVars) {
    const auto binnedDistance = binnedDistances.at(var.name());
    const auto covView = atlas::array::make_view<double, 2>(covariances[var.name()]);

    const size_t levels = var.getLevels();

    auto targetVarField = matchingOuterFspace.createField<double>(
                atlas::option::levels(levels) |
                atlas::option::name(var.name()) |
                atlas::option::halo(0));
    auto targetVarView = atlas::array::make_view<double, 2>(targetVarField);

    auto actualVarField = matchingOuterFspace.createField<double>(
                atlas::option::levels(levels) |
                atlas::option::name(var.name()) |
                atlas::option::halo(0));
    auto actualVarView = atlas::array::make_view<double, 2>(actualVarField);

    const auto lonlatView = atlas::array::make_view<double, 2>(innerFspace.lonlat());

    int matIndex = 0;
    for (atlas::idx_t jnode = 0; jnode < targetVarView.shape(0); jnode++) {
      // Extract jnode-th compressed row of the interpolationMatrix
      const size_t points = matOuter[jnode+1] - matOuter[jnode];
      std::vector<double> weights;
      std::vector<size_t> colIndex;
      weights.reserve(points);
      colIndex.reserve(points);
      while (matIndex < matOuter[jnode+1]) {
        weights.push_back(matData[matIndex]);
        colIndex.push_back(matInner[matIndex]);
        matIndex++;
      }

      // Extract matrix of distances D, typically of size 4 by 4.
      // D_ij is distance between points i and j on inner grid.
      std::vector<std::vector<double>> distances;
      for (size_t i = 0; i < points; i++) {
        const auto jnode1 = colIndex[i];
        const auto point1 = atlas::Point2(lonlatView(jnode1, atlas::LON),
                                          lonlatView(jnode1, atlas::LAT));
        std::vector<double> row;
        row.reserve(points);
        for (size_t j = 0; j < points; j++) {
          // Use symmetry to fill the distance matrix
          if (j < i) {
            // Lower part
            row.push_back(distances[j][i]);
          } else if (j == i) {
            // Diagonal
            row.push_back(0.0);
          } else {
            // Upper part
            const auto jnode2 = colIndex[j];
            const auto point2 = atlas::Point2(lonlatView(jnode2, atlas::LON),
                                              lonlatView(jnode2, atlas::LAT));
            const double dist = atlas::util::Earth().distance(point1, point2);
            row.push_back(dist);
          }
        }
        distances.push_back(row);
      }

      for (size_t lev = 0; lev < levels; lev++) {
        // Extract covariance matrix C associated to distance matrix
        std::vector<std::vector<double>> covs;
        std::vector<double> variances;
        variances.reserve(points);
        for (size_t i = 0; i < points; i++) {
          std::vector<double> row;
          row.reserve(points);
          for (size_t j = 0; j < points; j++) {
            if (j < i) {
              // Lower part
              row.push_back(covs[j][i]);
            } else {
              // Diagonal and upper part
              const double dist = distances[i][j];
              const double cov = lookUpCovarianceAtDistance(binnedDistance, covView, lev, dist);
              row.push_back(cov);
            }
          }
          covs.push_back(row);
          variances.push_back(row[i]);
        }

        // Compute target interpolated variance: v^T w
        // Where v is the vector of variances and w the interpolation weights
        const double targetInterpolatedVariance = std::inner_product(
                    weights.cbegin(), weights.cend(), variances.cbegin(), 0.0);
        targetVarView(jnode, lev) = targetInterpolatedVariance;

        // Compute variance of interpolated value: w^T C w
        // where C is the covariance matrix and w the interpolation weights
        std::vector<double> Cw;
        for (const auto & row : covs) {
          Cw.push_back(std::inner_product(
                    row.cbegin(), row.cend(), weights.cbegin(), 0.0));
        }
        const double varianceOfInterpolatedValue = std::inner_product(
                    weights.cbegin(), weights.cend(), Cw.cbegin(), 0.0);
        actualVarView(jnode, lev) = varianceOfInterpolatedValue;
      }
    }

    const double alpha = alphas[ivar];
    if (alpha > 0.0) {
      // Compute multiplicative correction
      auto field = matchingOuterFspace.createField<double>(
                    atlas::option::levels(levels) |
                    atlas::option::name(var.name()) |
                    atlas::option::halo(0));
      atlas::field::for_each_value(targetVarField,
                                   actualVarField,
                                   field,
                                   [&](const double target,
                                       const double actual,
                                       double & v){
               const double multTarget = alpha * target + (1 - alpha) * actual;
               v = std::sqrt(multTarget / actual);
                                   });
      multFset.add(field);
    }

    if (alpha != 1.0) {
      // Compute remaining additive correction
      auto field = matchingOuterFspace.createField<double>(
                    atlas::option::levels(levels) |
                    atlas::option::name(var.name()) |
                    atlas::option::halo(0));
      atlas::field::for_each_value(targetVarField,
                                   actualVarField,
                                   field,
                                   [&](const double target,
                                       const double actual,
                                       double & v){
               const double multTarget = alpha * target + (1 - alpha) * actual;
               v = std::sqrt(target - multTarget);
                                   });
      addFset.add(field);
    }

    ivar++;
  }

  oops::Log::trace() << "computeRescalingCoeffs about to exit..." << std::endl;
  return std::make_tuple(multFset, addFset);
}

// -----------------------------------------------------------------------------

atlas::FieldSet createRescalingCoeffs(
        const eckit::mpi::Comm & comm,
        const eckit::LocalConfiguration & conf,
        const oops::Variables & vars,
        const atlas::FunctionSpace & innerFspace,
        const atlas::FunctionSpace & outerFspace,
        const saber::interpolation::AtlasInterpWrapper & interp) {
  oops::Log::trace() << "createRescalingCoeffs starting" << std::endl;

  // Define list of variables to rescale
  oops::Variables activeVars(vars);
  if (conf.has("active variables")) {
    const auto & confVars = oops::Variables(conf.getStringVector("active variables"));
    ASSERT(confVars <= activeVars);
    activeVars.intersection(confVars);
  }

  // Check which fraction of the rescaling is applied multiplicatively
  const auto defaultAlpha = conf.getDouble("fraction of lost variance to restore multiplicatively",
                                           1.0);
  std::vector<double> alphas(activeVars.size(), defaultAlpha);

  // Optionally, the default value of alpha can be overwritten variable by variable
  const auto key = "list of fractions of lost variance to restore multiplicatively";
  if (conf.has(key)) {
    const auto alphaConf = conf.getSubConfiguration(key);
    int ivar = 0;
    for (const auto & var : activeVars) {
      if (alphaConf.has(var.name())) {
        alphas[ivar] = alphaConf.getDouble(var.name());
      }
      ivar++;
    }
  }

  // Read covariance profiles from files
  if (!conf.has("horizontal covariance profile file path")) {
    throw eckit::UserError("Please provide a covariance profile file path", Here());
  }
  const auto & filePath = conf.getString("horizontal covariance profile file path");
  const auto[distances, covariances] = readCovarianceProfiles(filePath, activeVars);

  // Create fields of rescaling factors
  auto[multFset, addFset] = computeRescalingCoeffs(activeVars, innerFspace, interp,
                                                   distances, covariances,
                                                   alphas);

  // Redistribute to outer functionSpace
  auto redistribute = [&](atlas::FieldSet & fset) {
    atlas::FieldSet outFset;
    for (const auto & field : fset) {
      outFset.add(outerFspace.createField<double>(atlas::option::name(field.name()) |
                                                  atlas::option::levels(field.levels())));
    }
    const auto & redistr = interp.getRedistribution();
    redistr.execute(fset, outFset);
    fset = outFset;
  };
  redistribute(multFset);
  redistribute(addFset);

  // Optionally, write additive rescaling coefficients
  if (conf.has("additive coefficients atlas file")) {
    const auto writeConf = conf.getSubConfiguration("additive coefficients atlas file");
    oops::Log::info() << "Writing additive coefficients to "
                      << writeConf.getString("filepath") << ".nc" << std::endl;
    util::writeFieldSet(comm, writeConf, addFset);
  }

  oops::Log::trace() << "createRescalingCoeffs about to exit..." << std::endl;
  return multFset;
}

// -----------------------------------------------------------------------------

atlas::FieldSet readRescalingCoeffs(const eckit::mpi::Comm & comm,
                                    const eckit::LocalConfiguration & conf,
                                    const oops::Variables & vars,
                                    const atlas::FunctionSpace & fspace) {
  oops::Log::trace() << "readRescalingCoeffs starting" << std::endl;

  // Define list of variables to rescale
  oops::Variables activeVars(vars);
  if (conf.has("active variables")) {
    const auto & confVars = oops::Variables(conf.getStringVector("active variables"));
    activeVars.intersection(confVars);
  }

  // Read input file
  atlas::FieldSet fset;
  util::readFieldSet(comm, fspace, activeVars, conf, fset);

  oops::Log::trace() << "readRescalingCoeffs about to exit..." << std::endl;
  return fset;
}

// -----------------------------------------------------------------------------
}  // namespace

namespace saber {
namespace interpolation {

// -----------------------------------------------------------------------------

Rescaling::Rescaling(const eckit::mpi::Comm & comm,
                     const eckit::LocalConfiguration & conf,
                     const oops::Variables & vars,
                     const atlas::FunctionSpace & innerFspace,
                     const atlas::FunctionSpace & outerFspace,
                     const saber::interpolation::AtlasInterpWrapper & interp) :
  multiplicativeRescalingCoeffs_(createRescalingCoeffs(comm, conf, vars,
                                         innerFspace, outerFspace, interp)) {
  oops::Log::trace() << classname() << "Rescaling starting" << std::endl;

  // Optionally, write rescaling fields to file
  if (conf.has("multiplicative coefficients atlas file")) {
    const auto writeConf = conf.getSubConfiguration("multiplicative coefficients atlas file");
    oops::Log::info() << "Writing multiplicative coefficients to "
                      << writeConf.getString("filepath") << ".nc" << std::endl;
    util::writeFieldSet(comm, writeConf, multiplicativeRescalingCoeffs_);
  }
  oops::Log::trace() << classname() << "Rescaling done" << std::endl;
}

// -----------------------------------------------------------------------------

Rescaling::Rescaling(const eckit::mpi::Comm & comm,
                     const eckit::LocalConfiguration & conf,
                     const oops::Variables & vars,
                     const atlas::FunctionSpace & outerFspace) :
    multiplicativeRescalingCoeffs_(readRescalingCoeffs(comm, conf, vars, outerFspace)) {
  oops::Log::trace() << classname() << "Rescaling done, from input file" << std::endl;
}

// -----------------------------------------------------------------------------

void Rescaling::execute(oops::FieldSet3D & fieldSet) const {
  oops::Log::trace() << classname() << "::execute starting" << std::endl;
  for (const auto & wField : multiplicativeRescalingCoeffs_) {
    const auto & field = fieldSet.fieldSet()[wField.name()];
    const auto & ghost = field.functionspace().ghost();
    atlas::field::for_each_value_masked(ghost,
                                        wField,
                                        field,
                                        [&](const double w, double & v){v *= w;});
    field.set_dirty();
  }
  oops::Log::trace() << classname() << "::execute done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace interpolation
}  // namespace saber
