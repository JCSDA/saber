/*
 * (C) Crown Copyright 2023-2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/generic/VertLoc.h"

#include <netcdf.h>

#include <iomanip>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/parallel/omp/omp.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/AtlasArrayUtil.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/oops/Utilities.h"

#define ERR(e) {throw eckit::Exception(nc_strerror(e), Here());}

namespace {

auto setInnerVars(const oops::Variables & outerVars,
                  const oops::Variables & activeVars,
                  const int nmods) {
  // Return variables which are a copy of outerVars, except for activeVars
  // which should have levels changed to nmods
  oops::Variables innerVars(outerVars);
  for (auto & var : activeVars) {
    innerVars[var.name()].setLevels(nmods);
  }
  return innerVars;
}

}  // namespace

namespace saber {
namespace vader {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<VertLoc> makerVertLoc_("mo_vertical_localization");

// -----------------------------------------------------------------------------

VertLoc::VertLoc(const oops::GeometryData & outerGeometryData,
                 const oops::Variables & outerVars,
                 const eckit::Configuration & covarConf,
                 const Parameters_ & params,
                 const oops::FieldSet3D & xb,
                 const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    innerGeometryData_(outerGeometryData),
    activeVars_(getActiveVars(params, outerVars)),
    nlevs_(activeVars_[0].getLevels()),
    nmods_(params.truncation.value()),
    innerVars_(setInnerVars(outerVars, activeVars_, nmods_)),
    ncfilepath_(params.VertLocParams.value().covStatFileName.value().get_value_or("")),
    locfilepath_(params.VertLocParams.value().locFileName.value()),
    locFieldName_(params.VertLocParams.value().locFieldName.value()),
    meanPressFieldName_(params.VertLocParams.value().meanPressFieldName.value().get_value_or("")),
    Umatrix_(Eigen::MatrixXd::Identity(nlevs_, nmods_))
{
  oops::Log::trace() << classname() << "::VertLoc starting" << std::endl;

  // Check all active variables have same number of levels
  for (const auto & var : activeVars_) {
    if (var.getLevels() != nlevs_) {
      oops::Log::error() << "Error    : Cannot deal with multiple vertical resolutions."
                         << std::endl;
      oops::Log::error() << "Error    : Vertical localization matrix has "
                         << nlevs_ << " levels, but variable " << var << " has "
                         << var.getLevels() << " levels." << std::endl;
      throw eckit::Exception("Inconsistent number of levels in VertLoc", Here());
    }
  }

  // Check range of nmods_
  if (nmods_ <= 1 || nmods_ > nlevs_) {
    oops::Log::error() << "Error    : Truncation should be between 1 and number"
                       << " of model levels (" << nlevs_ << ") but is "
                       << nmods_ << ". " << std::endl;
    throw eckit::BadParameter("Invalid value for number of vertical modes", Here());
  }

  // Read localization files and perform computations on root PE
  const std::size_t root(0);
  if (innerGeometryData_.comm().rank() == root) {
    // Read mean pressure values at alternate vertical stagger to localizaton matrix
    auto Fld1 = atlas::Field(meanPressFieldName_,
              atlas::array::make_datatype<double>(),
              atlas::array::make_shape(nlevs_+1));
    auto fv1 = atlas::array::make_view<double, 1>(Fld1);
    fv1.assign(0.0);

    VertLoc::readPressVec(ncfilepath_, meanPressFieldName_, fv1);

    std::vector<double> p_bar_mean;
    for (int i = 0; i < nlevs_+1; ++i) {
      p_bar_mean.push_back(fv1(i));
    }

    // Define an inner product W to weigh levels according to air mass
    std::vector<double> SqrtInnerProduct(nlevs_, 1.0);
    if (*std::max_element(std::begin(p_bar_mean),
                          std::end(p_bar_mean)) > 0.0) {
      std::vector<double> InnerProduct(nlevs_, 1.0);
      for (int i = 0; i < nlevs_; ++i) {
        InnerProduct[i] = p_bar_mean[i] - p_bar_mean[i+1];
      }

      const double SumInnerProduct = std::accumulate(InnerProduct.begin(),
                                                     InnerProduct.end(),
                                                     0.0);
      const double averageInnerProduct = SumInnerProduct / static_cast<double>(nlevs_);

      for (int i = 0; i < nlevs_; ++i) {
        InnerProduct[i] /= averageInnerProduct;
        SqrtInnerProduct[i] = std::sqrt(InnerProduct[i]);
      }
    } else {
      oops::Log::info() << "Info     : No valid mean pressure values at theta levels found"
                        << std::endl;
    }

    // Read localization matrix C
    auto Fld = atlas::Field(locFieldName_,
              atlas::array::make_datatype<double>(),
              atlas::array::make_shape(nlevs_, nlevs_));
    auto fv = atlas::array::make_view<double, 2>(Fld);
    fv.assign(0.0);

    VertLoc::readLocMat(locfilepath_, locFieldName_, fv);

    // Check localization has unit diagonal
    const bool enforceUnitDiagonal = !params.allowNonUnitDiagonal;
    const bool renormalizeDiagonal = params.renormalizeDiagonal;
    if (enforceUnitDiagonal) {
      constexpr double absoluteTolerance = 1.0e-6;
      for (int i = 0; i <nlevs_; ++i) {
        if (std::abs(fv(i, i) - 1.0) > absoluteTolerance) {
          std::stringstream message;
          message << "The vertical localization matrix prescribed has a"
                  << " non unit diagonal: " << i + 1<< "th element is "
                  << std::setprecision(15) << fv(i, i) << ". ";
          if (renormalizeDiagonal) {
            oops::Log::info() << "Info     : " << message.str() << std::endl;
            oops::Log::info() << "Info     : I now normalize by multiplying "
                              << "associated row and column. " << std::endl;
            const double normalization = std::sqrt(1.0 / fv(i, i));
            for (int j = 0; j < nlevs_; ++j) {
              fv(i, j) *= normalization;
              fv(j, i) *= normalization;
            }
          } else {
            oops::Log::error() << "Error    : " << message.str() << std::endl;
            oops::Log::error() << "Error    : If you want to use it as such, add key `reproduce "
                               << "bug non-unit diagonal: true` in the saber block configuration. "
                               << std::endl;
            oops::Log::error() << "Error    : If you want to perform an internal renormalization, "
                               << "add key `renormalize to unit diagonal: true` "
                               << "in the saber block configuration. "
                               << std::endl;
          throw eckit::UserError("Vertical localization matrix has non-unit diagonal in "
                                 + classname(), Here());
          }
        }
      }
    }

    Eigen::MatrixXd m = Eigen::MatrixXd::Zero(nlevs_, nlevs_);
    for (int i = 0; i < nlevs_; ++i) {
      for (int j = 0; j < nlevs_; ++j) {
        m(i, j) = fv(i, j);
      }
    }

    if (m.any() > 0.0) {
      // Apply inner product: C -> C' = W^{1/2} C W^{1/2}
      for (int i = 0; i < nlevs_; ++i) {
        for (int j = 0; j < nlevs_; ++j) {
          m(i, j) *= SqrtInnerProduct[i]*SqrtInnerProduct[j];
        }
      }

      // Do EOF decomposition of target matrix to get eigenvalues and
      // eigenvectors.
      //
      // C' == PDP^T with D diagonal and P unitary
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(m);

      // The eigenvalues are sorted in increasing order.
      // Reorder eigenvalues and eigenvectors to be in descending order
      Eigen::VectorXd evals(es.eigenvalues());
      Eigen::MatrixXd evecs(es.eigenvectors());
      for (int i = 0; i < nlevs_; ++i) {
        evals(i) = es.eigenvalues()(nlevs_-1-i);
        evecs.col(i) = es.eigenvectors().col(nlevs_-1-i);
      }

      // Calculate U := W^{-1/2} P D^{1/2} so that UU^T == C
      for (int j = 0; j < nmods_; ++j) {
        if (evals(j) > 0.0) {
          for (int i = 0; i < nlevs_; ++i) {
            evecs(i, j) *= std::sqrt(evals(j))/SqrtInnerProduct[i];
          }
        } else {
          for (int i = 0; i < nlevs_; ++i) {
            evecs(i, j) = 0.0;
          }
        }
      }

      // Total variance in weighted space
      const double totalVariance = std::accumulate(evals.begin(), evals.end(), 0.0);

      // Only keep first nmods_ modes into U
      double retainedVariance(0.0);
      for (int i = 0; i < nmods_; ++i) {
        Umatrix_.col(i) = evecs.col(i);
        retainedVariance += evals(i);
      }

      // Print percentage of variance retained
      const double percentageRetained = 100 * retainedVariance / totalVariance;
      oops::Log::info() << "Info     : Percentage of explained variance in "
                        << "pressure-weighted space is "
                        << std::setprecision(4) << percentageRetained << "%."
                        << std::endl;

      // Restore variance to original value at each level:
      double sumUmatrix(0.0);
      for (int i = 0; i < nlevs_; ++i) {
        if (m(i, i) > 0.0) {
          sumUmatrix = 0.0;
          for (int j = 0; j < nmods_; ++j) {
            sumUmatrix += Umatrix_(i, j)*Umatrix_(i, j);
          }
          for (int j = 0; j < nmods_; ++j) {
            Umatrix_(i, j) *= std::sqrt(fv(i, i))/std::sqrt(sumUmatrix);
          }
        } else {
          for (int j = 0; j < nmods_; ++j) {
            Umatrix_(i, j) = 0.0;
          }
        }
      }

      // Dump to file
      const auto & outputFileName = params.outputFileName.value();
      if (outputFileName != boost::none) {
        this->writeLocalization(outputFileName.value(),
                                SqrtInnerProduct,
                                fv,
                                Umatrix_);
      }

    } else {
      oops::Log::error() << "Error    : No valid localization matrix found in "
                         << locfilepath_
                         << std::endl;
      throw eckit::Exception("No valid localization matrix.", Here());
    }
  }
  innerGeometryData_.comm().broadcast(Umatrix_.data(), nmods_*nlevs_, root);

  oops::Log::trace() << classname() << "::VertLoc done" << std::endl;
}

// -----------------------------------------------------------------------------

VertLoc::~VertLoc() {
  oops::Log::trace() << classname() << "::~VertLoc starting" << std::endl;
  util::Timer timer(classname(), "~VertLoc");
  oops::Log::trace() << classname() << "::~VertLoc done" << std::endl;
}

// -----------------------------------------------------------------------------

void VertLoc::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  atlas::FieldSet fsetOut;

  // Passive variables
  for (const auto & var : fset.field_names()) {
    if (!activeVars_.has(var)) {
      fsetOut.add(fset[var]);
    }
  }

  // Active variables
  for (const auto & var : activeVars_) {
    if (fset[var].shape(1) != nmods_) {
      oops::Log::error() << "Error    : Field " << var << " has " << fset[var.name()].shape(1)
                         << ", expected " << nmods_ << ". " << std::endl;
      throw eckit::UserError("Wrong number of vertical levels in field " + var.name(), Here());
    }

    // Create new field with nlevs_ levels
    atlas::Field outField =
      innerGeometryData_.functionSpace().createField<double>
        (atlas::option::name(var.name()) |
         atlas::option::levels(nlevs_));
    auto outView = atlas::array::make_view<double, 2>(outField);
    outView.assign(0.0);

    // Apply U matrix
    auto inView = atlas::array::make_view<double, 2>(fset[var.name()]);  // nmods_ levels
    atlas_omp_parallel_for(atlas::idx_t jn = 0; jn < outField.shape(0); ++jn) {
      for (atlas::idx_t jl = 0; jl < nlevs_; ++jl) {
        for (atlas::idx_t jm = 0; jm < nmods_; ++jm) {
          outView(jn, jl) += Umatrix_(jl, jm) * inView(jn, jm);
        }
      }
    }

    fsetOut.add(outField);
  }

  fset.fieldSet() = fsetOut;

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void VertLoc::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  atlas::FieldSet fsetOut;

  // Passive variables
  for (const auto & var : fset.field_names()) {
    if (!activeVars_.has(var)) {
      fsetOut.add(fset[var]);
    }
  }

  // Active variables
  for (const auto & var : activeVars_) {
    if (fset[var.name()].shape(1) != nlevs_) {
      oops::Log::error() << "Error    : Field " << var << " has " << fset[var.name()].shape(1)
                         << ", expected " << nlevs_ << ". " << std::endl;
      throw eckit::UserError("Wrong number of vertical levels in field " + var.name(), Here());
    }

    // Create new field with nmods_ levels
    atlas::Field outField =
      innerGeometryData_.functionSpace().createField<double>
        (atlas::option::name(var.name()) |
         atlas::option::levels(nmods_));
    auto outView = atlas::array::make_view<double, 2>(outField);

    outView.assign(0.0);

    // Apply U^t
    auto inView = atlas::array::make_view<double, 2>(fset[var.name()]);  // nlevs_ levels

    for (atlas::idx_t jn = 0; jn < outField.shape(0); ++jn) {
      for (atlas::idx_t jl = 0; jl < nmods_; ++jl) {
        for (atlas::idx_t jm = 0; jm < nlevs_; ++jm) {
          outView(jn, jl) += Umatrix_(jm, jl) * inView(jn, jm);
        }
      }
    }


    fsetOut.add(outField);
  }

  fset.fieldSet() = fsetOut;

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void VertLoc::leftInverseMultiply(oops::FieldSet3D & fset) const {
  throw eckit::NotImplemented("leftInverseMultiply not implemented", Here());
}

// -----------------------------------------------------------------------------
void VertLoc::readLocMat(const std::string & filepath,
                         const std::string & fieldname,
                         atlas::array::ArrayView<double, 2> & fview) {
  // Setup
  std::vector<std::string> dimNames;
  std::vector<atlas::idx_t> dimSizes;
  std::vector<std::vector<std::string>> dimNamesForEveryVar;
  oops::Variables vars;
  std::vector<int> netcdfGeneralIDs;
  std::vector<int> netcdfDimIDs;
  std::vector<int> netcdfVarIDs;
  std::vector<std::vector<int>> netcdfDimVarIDs;
  eckit::LocalConfiguration netcdfMetaData;

  util::atlasArrayInquire(filepath,
                          dimNames,
                          dimSizes,
                          vars,
                          dimNamesForEveryVar,
                          netcdfMetaData,
                          netcdfGeneralIDs,
                          netcdfDimIDs,
                          netcdfVarIDs,
                          netcdfDimVarIDs);

  if (vars.has(fieldname)) {
    const std::size_t pos = vars.find(fieldname);
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
    if (dimFld != 2) throw eckit::Exception(fieldname + " is not a 2D field", Here());
    // check dimFldSizes[0] is nlevs_+1
    if (vecNumLevs != nlevs_) throw eckit::Exception(fieldname+" should have "
      + std::to_string(nlevs_) + " levels", Here());

    util::atlasArrayReadData(netcdfGeneralIDs,
                             netcdfVarID,
                             fview);

    int retval;
    if ((retval = nc_close(netcdfGeneralIDs[0]))) ERR(retval);
  }
}

// -----------------------------------------------------------------------------

void VertLoc::readPressVec(const std::string & filepath,
                           const std::string & fieldname,
                           atlas::array::ArrayView<double, 1> & fview) {
  // Setup
  std::vector<std::string> dimNames;
  std::vector<atlas::idx_t> dimSizes;
  oops::Variables vars;
  std::vector<std::vector<std::string>> dimNamesForEveryVar;
  eckit::LocalConfiguration netcdfMetaData;
  std::vector<int> netcdfGeneralIDs;
  std::vector<int> netcdfDimIDs;
  std::vector<int> netcdfVarIDs;
  std::vector<std::vector<int>> netcdfDimVarIDs;

  std::string nullstr = "";
  if (ncfilepath_.compare(nullstr) > 0) {
    util::atlasArrayInquire(ncfilepath_,
                            dimNames,
                            dimSizes,
                            vars,
                            dimNamesForEveryVar,
                            netcdfMetaData,
                            netcdfGeneralIDs,
                            netcdfDimIDs,
                            netcdfVarIDs,
                            netcdfDimVarIDs);

    if (meanPressFieldName_.compare(nullstr) > 0) {
      if (vars.has(meanPressFieldName_)) {
        const std::size_t pos = vars.find(meanPressFieldName_);
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
        if (dimFld > 1) throw eckit::Exception(meanPressFieldName_ + " is not a 1D field",
                                               Here());
        // check dimFldSizes[0] is nlevs_+1
        if (vecNumLevs != nlevs_+1) throw eckit::Exception(meanPressFieldName_ + " should have " +
        std::to_string(nlevs_+1) + " levels", Here());

        util::atlasArrayReadData(netcdfGeneralIDs,
                                 netcdfVarID,
                                 fview);
      }
    }
    int retval;
    if ((retval = nc_close(netcdfGeneralIDs[0]))) ERR(retval);
  }
}

// -----------------------------------------------------------------------------

void VertLoc::writeLocalization(
        const std::string & filepath,
        const std::vector<double> & sqrtInnerProduct,
        const atlas::array::ArrayView<const double, 2> & targetLocalizationView,
        const Eigen::MatrixXd & localizationSquareRoot) {
  // This function is meant to run on a single PE
  oops::Log::trace() << classname() << "::writeLocalization starting" << std::endl;

  // Compute low rank localization
  const auto lowRankLoc = localizationSquareRoot * localizationSquareRoot.transpose();

  // Convert all data holders to atlas fields
  atlas::FieldSet fset;

  auto field0 = atlas::Field(
              "air_mass_weights",
              atlas::array::make_datatype<double>(),
              atlas::array::make_shape(nlevs_));
  auto view0 = atlas::array::make_view<double, 1>(field0);
  for (int i = 0; i < nlevs_; ++i) {
    view0(i) = sqrtInnerProduct[i];
  }
  fset.add(field0);

  auto field1 = atlas::Field("target_localization",
                            atlas::array::make_datatype<double>(),
                            atlas::array::make_shape(nlevs_, nlevs_));
  auto view1 = atlas::array::make_view<double, 2>(field1);
  view1.assign(targetLocalizationView);
  fset.add(field1);

  auto field2 = atlas::Field(
              "localization_square_root",
              atlas::array::make_datatype<double>(),
              atlas::array::make_shape(nlevs_, nmods_));
  auto view2 = atlas::array::make_view<double, 2>(field2);
  for (int i = 0; i < nlevs_; ++i) {
    for (int j = 0; j < nmods_; ++j) {
      view2(i, j) = localizationSquareRoot(i, j);
    }
  }
  fset.add(field2);

  auto field3 = atlas::Field(
              "low_rank_localization",
              atlas::array::make_datatype<double>(),
              atlas::array::make_shape(nlevs_, nlevs_));
  auto view3 = atlas::array::make_view<double, 2>(field3);
  for (int i = 0; i < nlevs_; ++i) {
    for (int j = 0; j < nlevs_; ++j) {
      view3(i, j) = lowRankLoc(i, j);
    }
  }
  fset.add(field3);

  // Define variables and dimensions
  const std::vector<std::string> dimNames{"nz", "nmods"};
  const std::vector<atlas::idx_t> dimSizes{nlevs_, nmods_};

  // variables constructed from vector of strings to initialize
  // empty metadata. This may be changed in future.
  const oops::Variables vars(std::vector<std::string>{"air_mass_weights",
                                                      "target_localization",
                                                      "low_rank_localization",
                                                      "localization_square_root"});
  const std::vector<std::vector<std::string>>
    dimNamesForEveryVar{{"nz"},
                        {"nz", "nz"},
                        {"nz", "nz"},
                        {"nz", "nmods"}};
  const std::vector<std::vector<atlas::idx_t>>
    dimSizesForEveryVar{{nlevs_},
                        {nlevs_, nlevs_},
                        {nlevs_, nlevs_},
                        {nlevs_, nmods_}};

  eckit::LocalConfiguration netcdfMetaData;
  for (const oops::Variable & var : vars) {
    util::setAttribute<std::string>(
      netcdfMetaData, var.name(), "statistics_type", "string", "vertical localization");
  }

  // Write Header
  std::vector<int> netcdfGeneralIDs;
  std::vector<int> netcdfDimIDs;
  std::vector<int> netcdfVarIDs;
  std::vector<std::vector<int>> netcdfDimVarIDs;

  util::atlasArrayWriteHeader(filepath,
                              dimNames,
                              dimSizes,
                              vars,
                              dimNamesForEveryVar,
                              netcdfMetaData,
                              netcdfGeneralIDs,
                              netcdfDimIDs,
                              netcdfVarIDs,
                              netcdfDimVarIDs);

  // Write Data
  // Needs new views with datatype `const double`
  int t = 0;
  for (const auto & var : vars) {
    const auto field = fset.field(var.name());
    const size_t rank = field.shape().size();
    if (rank == 1) {
      auto view = atlas::array::make_view<const double, 1>(field);
      util::atlasArrayWriteData(netcdfGeneralIDs, netcdfVarIDs[t], view);
    } else if (rank == 2) {
      auto view = atlas::array::make_view<const double, 2>(field);
      util::atlasArrayWriteData(netcdfGeneralIDs, netcdfVarIDs[t], view);
    }
    t++;
  }

  oops::Log::info() << "Info     : Localization data written to file "
                    << filepath << std::endl;

  int retval;
  if ((retval = nc_close(netcdfGeneralIDs[0]))) throw eckit::Exception("NetCDF closing error",
                                                                       Here());

  oops::Log::trace() << classname() << "::writeLocalization done" << std::endl;
}
// -----------------------------------------------------------------------------

void VertLoc::print(std::ostream & os) const {
  os << classname();
}

}  // namespace vader
}  // namespace saber
