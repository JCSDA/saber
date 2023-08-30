/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/vader/VertLoc.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/AtlasArrayUtil.h"
#include "oops/util/Timer.h"

#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/vader/CovarianceStatisticsUtils.h"

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<VertLoc> makerVertLoc_("mo_vertical_localisation");

// -----------------------------------------------------------------------------

VertLoc::VertLoc(const oops::GeometryData & outerGeometryData,
                                   const oops::Variables & outerVars,
                                   const eckit::Configuration & covarConf,
                                   const Parameters_ & params,
                                   const oops::FieldSet3D & xb,
                                   const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params),
    innerGeometryData_(outerGeometryData), innerVars_(outerVars),
    activeVars_(params.activeVars.value().get_value_or(outerVars)),
    nlevs(params.nlevels.value()),
    nmods(params.truncation.value()),
    ncfilepath_(params.VertLocParams.value().covStatFileName.value().get_value_or("")),
    locfilepath_(params.VertLocParams.value().locFileName.value()),
    locFieldName_(params.VertLocParams.value().locFieldName.value()),
    meanPressFieldName_(params.VertLocParams.value().meanPressFieldName.value().get_value_or("")),
    Umatrix_(Eigen::MatrixXd::Identity(nlevs, nmods))
{
  oops::Log::trace() << classname() << "::VertLoc starting" << std::endl;

  // read from file on root PE if required.
  std::size_t root(0);
  if (innerGeometryData_.comm().rank() == root) {
    auto Fld1 = atlas::Field(meanPressFieldName_,
              atlas::array::make_datatype<double>(),
              atlas::array::make_shape(nlevs+1));
    auto fv1 = atlas::array::make_view<double, 1>(Fld1);
    fv1.assign(0.0);

    // read mean pressure values at theta levels
    VertLoc::readPressVec(ncfilepath_, meanPressFieldName_, fv1);

    std::vector<double> ptheta_bar_mean;
    for (int i = 0; i < nlevs+1; ++i) {
      ptheta_bar_mean.push_back(fv1(i));
    }

    std::vector<double> SqrtInnerProduct(nlevs, 1.0);
    if (*std::max_element(std::begin(ptheta_bar_mean),
                          std::end(ptheta_bar_mean)) > 0.0) {
      std::vector<double> InnerProduct(nlevs, 1.0);
      for (int i = 0; i < nlevs; ++i) {
        InnerProduct[i] = ptheta_bar_mean[i] - ptheta_bar_mean[i+1];
      }

      double SumInnerProduct(0.0);
      for (int i = 0; i < nlevs; ++i) {
        SumInnerProduct += InnerProduct[i];
      }

      // normalise inner product by the average delta pressure
      // before taking the square root
      for (int i = 0; i < nlevs; ++i) {
        InnerProduct[i] *= static_cast<double>(nlevs);
        InnerProduct[i] /= SumInnerProduct;
        SqrtInnerProduct[i] = std::sqrt(InnerProduct[i]);
      }
    } else {
      oops::Log::info() << "No valid mean pressure values at theta levels found"
                        << std::endl;
    }

    auto Fld = atlas::Field(locFieldName_,
              atlas::array::make_datatype<double>(),
              atlas::array::make_shape(nlevs, nlevs));
    auto fv = atlas::array::make_view<double, 2>(Fld);
    fv.assign(0.0);

    // read localisation matrix
    VertLoc::readLocMat(locfilepath_, locFieldName_, fv);

    Eigen::MatrixXd m = Eigen::MatrixXd::Zero(nlevs, nlevs);

    for (int i = 0; i < nlevs; ++i) {
      for (int j = 0; j < nlevs; ++j) {
        m(i, j) = fv(i, j);
      }
    }

    if (m.any() > 0.0) {
      for (int i = 0; i < nlevs; ++i) {
        for (int j = 0; j < nlevs; ++j) {
          m(i, j) *= SqrtInnerProduct[i]*SqrtInnerProduct[j];
        }
      }

      // Do EOF decomposition of target matrix to get eigenvalues and
      // eigenvectors.
      //
      // C --> PCP

      // The eigenvalues are sorted in increasing order
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(m);

      // reorder eigenvalues and eigenvectors to be in descending order
      Eigen::VectorXd evals(es.eigenvalues());
      Eigen::MatrixXd evecs(es.eigenvectors());
      for (int i = 0; i < nlevs; ++i) {
        evals(i) = es.eigenvalues()(nlevs-1-i);
        evecs.col(i) = es.eigenvectors().col(nlevs-1-i);
      }

      // calculates D := P^{-1} E \lambda^{1/2} so that DD^T = C
      for (int j = 0; j < nlevs; ++j) {
        if (evals(j) > 0.0) {
          for (int i = 0; i < nlevs; ++i) {
            evecs(i, j) *= std::sqrt(evals(j))/SqrtInnerProduct[i];
          }
        } else {
          for (int i = 0; i < nlevs; ++i) {
            evecs(i, j) = 0.0;
          }
        }
      }

      for (int i = 0; i < nmods; ++i) {
        Umatrix_.col(i) = evecs.col(i);
      }

      // Restore variance to target value at each level:
      double sumUmatrix(0.0);
      for (int i = 0; i < nlevs; ++i) {
        if (m(i, i) > 0.0) {
          sumUmatrix = 0.0;
          for (int j = 0; j < nmods; ++j) {
            sumUmatrix += Umatrix_(i, j)*Umatrix_(i, j);
          }
          for (int j = 0; j < nmods; ++j) {
            Umatrix_(i, j) *= std::sqrt(fv(i, i))/std::sqrt(sumUmatrix);
          }
        } else {
          for (int j = 0; j < nmods; ++j) {
            Umatrix_(i, j) = 0.0;
          }
        }
      }
    } else {
      oops::Log::info() << "No valid localisation matrix found"
                          << std::endl;
    }
  }
  innerGeometryData_.comm().broadcast(Umatrix_.data(), nmods*nlevs, root);

  oops::Log::trace() << classname() << "::VertLoc done" << std::endl;
}

// -----------------------------------------------------------------------------

VertLoc::~VertLoc() {
  oops::Log::trace() << classname() << "::~VertLoc starting" << std::endl;
  util::Timer timer(classname(), "~VertLoc");
  oops::Log::trace() << classname() << "::~VertLoc done" << std::endl;
}

// -----------------------------------------------------------------------------

void VertLoc::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  oops::Log::trace() << "active variables: " << activeVars_.variables() << std::endl;

  if (nmods != nlevs)
    ABORT("Number of vertical modes is here supposed to be equal to number of vertical levels");

  for (const auto & var : activeVars_.variables()) {
    auto fcsView = atlas::array::make_view<double, 2>(fset[var]);  // (nlocs, nmods)
    std::vector<double> kRowByCol(nmods);
    for (atlas::idx_t jn = 0; jn < fset[var].shape(0); ++jn) {
      for (atlas::idx_t jl = 0; jl < nlevs; ++jl) {
        kRowByCol[jl] = 0.0;
        for (atlas::idx_t jm = 0; jm < nmods; ++jm) {
          kRowByCol[jl] += Umatrix_(jl, jm) * fcsView(jn, jm);
        }
      }
      for (atlas::idx_t jl = 0; jl < nlevs; ++jl) {
        // this is the (nlocs, nlevs) localised field
        // for now output same as input (must be nmods = nlevs)
        fcsView(jn, jl) =  kRowByCol[jl];
      }
    }
  }

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void VertLoc::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  if (nmods != nlevs)
    ABORT("Number of vertical modes is here supposed to be equal to number of vertical levels");

  for (const auto & var : activeVars_.variables()) {
    auto fcsHatView =
      atlas::array::make_view<double, 2>(fset[var]);  // (nlocs, nlevs)
    std::vector<double> kColByRow(nmods);
    for (atlas::idx_t jn = 0; jn < fset[var].shape(0); ++jn) {
      for (atlas::idx_t jl = 0; jl < nmods; ++jl) {
        kColByRow[jl] = 0.0;
        for (atlas::idx_t jm = 0; jm < nlevs; ++jm) {
          kColByRow[jl] += Umatrix_(jm, jl) * fcsHatView(jn, jm);
        }
      }
      for (atlas::idx_t jl = 0; jl < nmods; ++jl) {
        // for now output same as input (must be nmods = nlevs)
        // note here we set = and not += given input field same as output
        fcsHatView(jn, jl) = kColByRow[jl];
      }
    }
  }

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void VertLoc::leftInverseMultiply(atlas::FieldSet & fset) const {
  throw eckit::NotImplemented(
          "leftInverseMultiply not implemented", Here());
}

// -----------------------------------------------------------------------------
void VertLoc::readLocMat(const std::string filepath, const std::string fieldname,
                         atlas::array::ArrayView<double, 2> & fview) {
  // Setup
  std::vector<std::string> dimNames;
  std::vector<atlas::idx_t> dimSizes;
  std::vector<std::vector<std::string>> dimNamesForEveryVar;
  std::vector<std::string> variableNames;
  std::vector<int> netcdfGeneralIDs;
  std::vector<int> netcdfDimIDs;
  std::vector<int> netcdfVarIDs;
  std::vector<std::vector<int>> netcdfDimVarIDs;

  util::atlasArrayInquire(filepath,
                          dimNames,
                          dimSizes,
                          variableNames,
                          dimNamesForEveryVar,
                          netcdfGeneralIDs,
                          netcdfDimIDs,
                          netcdfVarIDs,
                          netcdfDimVarIDs);

  auto ind = std::find(variableNames.begin(), variableNames.end(), fieldname);
  if (ind != variableNames.end()) {
    auto pos = ind - variableNames.begin();
    auto netcdfVarID = netcdfVarIDs[pos];

    int dimFld = netcdfDimVarIDs[pos].size();
    // remove dummy dimensions
    for (auto netcdfDimVarID : netcdfDimVarIDs[pos]) {
      if (dimSizes[netcdfDimVarID] == 1) dimFld -= 1;
    }

    std::vector<int> dimFldSizes;
    for (int i = 0; i < dimFld; ++i) dimFldSizes.push_back(dimSizes[netcdfDimVarIDs[pos][i]]);
    auto vecNumLevs = *std::max_element(std::begin(dimFldSizes), std::end(dimFldSizes));

    // current code assumes field in netcdf file to be 2D
    if (dimFld != 2) ABORT(fieldname+" is not a 2D field");
    // check dimFldSizes[0] is nlevs+1
    if (vecNumLevs != nlevs) ABORT(fieldname+" should have "+std::to_string(nlevs)+" levels");

    util::atlasArrayReadData(netcdfGeneralIDs,
                             dimFldSizes,
                             netcdfVarID,
                             fview);
    int retval;
    if ((retval = nc_close(netcdfGeneralIDs[0]))) ABORT(nc_strerror(retval));
  }
}

// -----------------------------------------------------------------------------

void VertLoc::readPressVec(const std::string filepath, const std::string fieldname,
                         atlas::array::ArrayView<double, 1> & fview) {
  // Setup
  std::vector<std::string> dimNames;
  std::vector<atlas::idx_t> dimSizes;
  std::vector<std::vector<std::string>> dimNamesForEveryVar;
  std::vector<std::string> variableNames;
  std::vector<int> netcdfGeneralIDs;
  std::vector<int> netcdfDimIDs;
  std::vector<int> netcdfVarIDs;
  std::vector<std::vector<int>> netcdfDimVarIDs;

  std::string nullstr = "";
  if (ncfilepath_.compare(nullstr) > 0) {
    util::atlasArrayInquire(ncfilepath_,
                            dimNames,
                            dimSizes,
                            variableNames,
                            dimNamesForEveryVar,
                            netcdfGeneralIDs,
                            netcdfDimIDs,
                            netcdfVarIDs,
                            netcdfDimVarIDs);

    if (meanPressFieldName_.compare(nullstr) > 0) {
      auto ind = std::find(variableNames.begin(), variableNames.end(), meanPressFieldName_);
      if (ind != variableNames.end()) {
        auto pos = ind - variableNames.begin();
        auto netcdfVarID = netcdfVarIDs[pos];
        int dimFld = netcdfDimVarIDs[pos].size();
        // remove dummy dimensions
        for (auto netcdfDimVarID : netcdfDimVarIDs[pos]) {
          if (dimSizes[netcdfDimVarID] == 1) dimFld -= 1;
        }

        std::vector<int> dimFldSizes;
        for (int i = 0; i < dimFld; ++i) dimFldSizes.push_back(dimSizes[netcdfDimVarIDs[pos][i]]);
        auto vecNumLevs = *std::max_element(std::begin(dimFldSizes), std::end(dimFldSizes));

        // current code assumes field in netcdf file to be 1D
        if (dimFld > 1) ABORT(meanPressFieldName_+" is not a 1D field");
        // check dimFldSizes[0] is nlevs+1
        if (vecNumLevs != nlevs+1) ABORT(meanPressFieldName_+" should have "+
                                         std::to_string(nlevs+1)+" levels");

        util::atlasArrayReadData(netcdfGeneralIDs,
                                 dimFldSizes,
                                 netcdfVarID,
                                 fview);
        int retval;
        if ((retval = nc_close(netcdfGeneralIDs[0]))) ABORT(nc_strerror(retval));
      }
    }
  }
}

// -----------------------------------------------------------------------------

void VertLoc::print(std::ostream & os) const {
  os << classname();
}

}  // namespace vader
}  // namespace saber
