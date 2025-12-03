/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/BifourierCovariance.h"

#include <Eigen/Dense>
#include <netcdf.h>

#include <algorithm>

#include "oops/generic/gc99.h"

#include "saber/bifourier/BifourierUtilities.h"

#define ERR(e, msg) {std::string s(nc_strerror(e)); \
  throw eckit::Exception(s + " : " + msg, Here());}

using atlas::array::make_datatype;
using atlas::array::make_indexview;
using atlas::array::make_shape;
using atlas::array::make_view;

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

static SaberCentralBlockMaker<BifourierCovariance> makerBifourierCovariance_("BifourierCovariance");

// -----------------------------------------------------------------------------

BifourierCovariance::BifourierCovariance(const oops::GeometryData & gdata,
                                         const oops::Variables & activeVars,
                                         const eckit::Configuration & covarConf,
                                         const Parameters_ & params,
                                         const oops::FieldSet3D & xb,
                                         const oops::FieldSet3D & fg)
  : SaberCentralBlockBase(params, xb.validTime()),
    comm_(gdata.comm()),
    activeVars_(activeVars),
    params_(params),
    Lf_(params_.calibration.value() != boost::none ?
      params_.calibration.value()->filteringScale.value() : 0),
    trans_(transStore_.retrieveTransform(gdata))
{
  oops::Log::trace() << classname() << "::BifourierCovariance starting" << std::endl;
  oops::Log::trace() << classname() << "::BifourierCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

BifourierCovariance::~BifourierCovariance() {
  oops::Log::trace() << classname() << "::~BifourierCovariance starting" << std::endl;
  oops::Log::trace() << classname() << "::~BifourierCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierCovariance::randomize(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;

  // Create random spectral vector
  trans_->createRandomFieldSet(fset.fieldSet(), activeVars_);

  // Define control vector
  atlas::Field cv("cv", make_datatype<double>(), make_shape(ctlVecSize()));

  // No offset
  const size_t offset = 0;

  // Convert spectral FieldSet to control vector
  trans_->fset2cv(fset.fieldSet(), cv, activeVars_, offset);

  // Square-root multiply
  multiplySqrt(cv, fset, offset);

  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierCovariance::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  // Define control vector
  atlas::Field cv("cv", make_datatype<double>(), make_shape(ctlVecSize()));

  // No offset
  const size_t offset = 0;

  // Square-root adjoint multiply
  multiplySqrtAD(fset, cv, offset);

  // Square-root multiply
  multiplySqrt(cv, fset, offset);

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierCovariance::multiplySqrt(const atlas::Field & cv,
                                       oops::FieldSet3D & fset,
                                       const size_t & offset) const {
  oops::Log::trace() << classname() << "::multiplySqrt starting" << std::endl;

  // Convert control vector to spectral FieldSet
  trans_->cv2fset(cv, fset.fieldSet(), activeVars_, offset);

  for (const auto & var : activeVars_) {
    // Get number of levels
    const size_t nz = var.getLevels();

    // Copy field
    auto copyField = fset[var.name()].clone();

    // Get field views
    auto view = getView2D(var, fset);
    auto copyView = make_view<double, 2>(copyField);

    // Get correlation square-root view
    const auto corSqrtView = getView3D("corSqrt", var, data_);

    // Get standard-deviation view
    const auto stdDevView = getViewProfile("stdDev", var, data_);

    // Set output field to zero
    view.assign(0.0);

    // Apply correlation square-root
    for (size_t js = 0; js < trans_->ns(); ++js) {
      const size_t jw = trans_->jw(js);
      for (size_t jzI = 0; jzI < nz; ++jzI) {
        for (size_t jzJ = 0; jzJ < nz; ++jzJ) {
          view(js, jzI) += corSqrtView(jw, jzI, jzJ)*copyView(js, jzJ);
        }
      }
    }

    if (!params_.correlation.value()) {
      // Apply standard-deviation
      for (size_t js = 0; js < trans_->ns(); ++js) {
        for (size_t jz = 0; jz < nz; ++jz) {
          view(js, jz) *= stdDevView(jz)*params_.inflation.value();
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::multiplySqrt done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierCovariance::multiplySqrtAD(const oops::FieldSet3D & fset,
                                         atlas::Field & cv,
                                         const size_t & offset) const {
  oops::Log::trace() << classname() << "::multiplySqrtAD starting" << std::endl;

  // Temporary FieldSet3D
  oops::FieldSet3D fsetTmp(fset);

  // Copy FieldSet
  trans_->copyFieldSet(fset.fieldSet(), fsetTmp.fieldSet(), activeVars_);

  for (const auto & var : activeVars_) {
    // Get number of levels
    const size_t nz = var.getLevels();

    // Copy field
    auto copyField = fsetTmp[var.name()].clone();

    // Get field views
    auto view = getView2D(var, fsetTmp);
    auto copyView = make_view<double, 2>(copyField);

    // Get correlation square-root view
    const auto corSqrtView = getView3D("corSqrt", var, data_);

    // Get standard-deviation view
    const auto stdDevView = getViewProfile("stdDev", var, data_);

    // Set output field to zero
    view.assign(0.0);

    if (!params_.correlation.value()) {
      // Apply standard-deviation
      for (size_t js = 0; js < trans_->ns(); ++js) {
        for (size_t jzJ = 0; jzJ < nz; ++jzJ) {
          copyView(js, jzJ) *= stdDevView(jzJ)*params_.inflation.value();
        }
      }
    }

    // Apply covariance square-root
    for (size_t js = 0; js < trans_->ns(); ++js) {
      const size_t jw = trans_->jw(js);
      for (size_t jzI = 0; jzI < nz; ++jzI) {
        for (size_t jzJ = 0; jzJ < nz; ++jzJ) {
          view(js, jzI) += corSqrtView(jw, jzJ, jzI)*copyView(js, jzJ);
        }
      }
    }
  }

  // Convert spectral FieldSet to control vector
  trans_->fset2cv(fsetTmp.fieldSet(), cv, activeVars_, offset);

  oops::Log::trace() << classname() << "::multiplySqrtAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierCovariance::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;

  for (const auto & var : activeVars_) {
    // Create correlation square-root field
    createField3D("corSqrt", trans_->nw(), var, data_);

    // Create standard-deviation field
    createFieldProfile("stdDev", var, data_);
  }

  // NetCDF file path
  const std::string ncFilePath = params_.read.value()->inputFile.value();

  // NetCDF IDs
  int ncId, retval, nwId, nzIId, nzJId, corSqrtId, stdDevId;
  size_t nwGlbFromFile, nzIFromFile, nzJFromFile;

  if (comm_.rank() == 0) {
    // Open NetCDF file
    if ((retval = nc_open(ncFilePath.c_str(), NC_NOWRITE, &ncId))) ERR(retval, ncFilePath);
  }

  for (const auto & var : activeVars_) {
    // Get number of levels
    const size_t nz = var.getLevels();

    // Define global vectors
    std::vector<double> corSqrtVecGlb;
    std::vector<double> stdDevVec(nz);

    if (comm_.rank() == 0) {
      // Check dimensions
      const std::string nzIName = "nzI_" + var.name();
      const std::string nzJName = "nzJ_" + var.name();
      if ((retval = nc_inq_dimid(ncId, "nwGlb", &nwId))) ERR(retval, "nwGlb");
      if ((retval = nc_inq_dimid(ncId, nzIName.c_str(), &nzIId))) ERR(retval, nzIName);
      if ((retval = nc_inq_dimid(ncId, nzJName.c_str(), &nzJId))) ERR(retval, nzJName);
      if ((retval = nc_inq_dimlen(ncId, nwId, &nwGlbFromFile))) ERR(retval, "nwGlb");
      if ((retval = nc_inq_dimlen(ncId, nzIId, &nzIFromFile))) ERR(retval, nzIName);
      if ((retval = nc_inq_dimlen(ncId, nzJId, &nzJFromFile))) ERR(retval, nzJName);
      ASSERT(nwGlbFromFile == trans_->nwGlb());
      ASSERT(nzIFromFile == nz);
      ASSERT(nzJFromFile == nz);

      // Get correlation square-root field name
      const std::string corSqrtFieldName = fieldName("corSqrt", var);

      // Get standard-deviation field name
      const std::string stdDevFieldName = fieldName("stdDev", var);

      // Get variables ID
      if ((retval = nc_inq_varid(ncId, corSqrtFieldName.c_str(), &corSqrtId)))
        ERR(retval, corSqrtFieldName);
      if ((retval = nc_inq_varid(ncId, stdDevFieldName.c_str(), &stdDevId)))
        ERR(retval, stdDevFieldName);

      // Allocate global correlation square-root vector
      corSqrtVecGlb.resize(trans_->nwGlb()*nz*nz);

      // Read data
      if ((retval = nc_get_var_double(ncId, corSqrtId, corSqrtVecGlb.data())))
        ERR(retval, corSqrtFieldName);
      if ((retval = nc_get_var_double(ncId, stdDevId, stdDevVec.data())))
        ERR(retval, stdDevFieldName);
    }

    // Get correlation square-root field
    auto corSqrtField = getField("corSqrt", var, data_);

    // Scatter correlation square-root vector
    trans_->scatterCov(corSqrtVecGlb, corSqrtField);

    // Broadcast standard-deviation vector
    comm_.broadcast(stdDevVec.begin(), stdDevVec.end(), 0);

    // Get standard-deviation view
    auto stdDevView = getViewProfile("stdDev", var, data_);

    // Deserialize standard-deviation vector
    for (size_t jzJ = 0; jzJ < nz; ++jzJ) {
      stdDevView(jzJ) = stdDevVec[jzJ];
    }
  }

  if (comm_.rank() == 0) {
    // Close file
    if ((retval = nc_close(ncId))) ERR(retval, ncFilePath);
  }

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierCovariance::directCalibration(const oops::FieldSets & fsetEns) {
  oops::Log::trace() << classname() << "::directCalibration starting" << std::endl;

  // Check ensemble size
  const size_t ne = fsetEns.ens_size();

  if (ne == 0) {
    // User-defined vertical profile for each variable

    // Index fields views
    const atlas::functionspace::StructuredColumns fs(trans_->geometryData().functionSpace());
    const auto indexIView = make_indexview<int, 1>(fs.index_i());
    const auto indexJView = make_indexview<int, 1>(fs.index_j());

    // Define fieldsets
    atlas::FieldSet horCorGpFset;
    atlas::FieldSet horCorSpFset;

    for (const auto & var : activeVars_) {
      // Get number of levels
      const size_t nz = var.getLevels();

      // Get horizontal length-scale profile
      std::vector<double> Lh;
      for (const auto & profile : params_.calibration.value()->profiles.value()) {
        if (profile.variable.value() == var.name()) {
          // Check profile size
          ASSERT(profile.Lh.value().size() == nz);

          // Allocate profiles
          ASSERT(Lh.size() == 0);
          Lh.resize(nz);

          // Copy horizontal length-scale profile
          Lh = profile.Lh.value();
        }
      }
      ASSERT(Lh.size() == nz);

      // Create horizontal grid-point correlation field
      auto horCorField = fs.createField<double>(atlas::option::name(var.name())
        | atlas::option::levels(nz));

      // Get horizontal grid-point correlation view
      auto horCorView = make_view<double, 2>(horCorField);

      // Compute horizontal grid-point correlation
      for (int jnode = 0; jnode < fs.size(); ++jnode) {
        const double distI = indexIView(jnode) < static_cast<int>(trans_->nx()/2) ?
          static_cast<double>(indexIView(jnode))*trans_->dx() :
          static_cast<double>(trans_->nx()-indexIView(jnode))*trans_->dx();
        const double distJ = indexJView(jnode) < static_cast<int>(trans_->ny()/2) ?
          static_cast<double>(indexJView(jnode))*trans_->dy() :
          static_cast<double>(trans_->ny()-indexJView(jnode))*trans_->dy();
        const double dist = std::sqrt(distI*distI+distJ*distJ);
        for (size_t jzI = 0; jzI < nz; ++jzI) {
          const double normDist = dist/Lh[jzI];
          horCorView(jnode, jzI) = oops::gc99(normDist);
        }
      }

      // Add field
      horCorGpFset.add(horCorField);
    }

    // Direct spectral transform of the horizontal grid-point correlation
    trans_->gp2sp(horCorGpFset, horCorSpFset, activeVars_);

    for (const auto & var : activeVars_) {
      // Get number of levels
      const size_t nz = var.getLevels();

      // Create correlation square-root
      atlas::Field corSqrtField("corSqrt", make_datatype<double>(),
        make_shape(trans_->nw(), nz, nz));

      // Get correlation square-root view
      auto corSqrtView = make_view<double, 3>(corSqrtField);

      // Get horizontal length-scale profile and vertical length-scale
      double Lv = 0.0;
      for (const auto & profile : params_.calibration.value()->profiles.value()) {
        if (profile.variable.value() == var.name()) {
          // Copy vertical length-scale
          Lv = profile.Lv.value();
        }
      }
      ASSERT(Lv > 0.0);

      // Compute vertical correlation matrix
      Eigen::MatrixXd vertCor(nz, nz);
      for (size_t jzI = 0; jzI < nz; ++jzI) {
        for (size_t jzJ = 0; jzJ < nz; ++jzJ) {
          const double normDist = std::abs(static_cast<double>(jzI)-static_cast<double>(jzJ))/Lv;
          vertCor(jzI, jzJ) = oops::gc99(normDist);
        }
      }

      // Compute eigendecomposition of the vertical correlation matrix
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
      es.compute(vertCor);

      // Get horizontal correlation in spectral space
      const auto horCorSpView = make_view<double, 2>(horCorSpFset[var.name()]);

      // Create horizontal spectral variance field
      atlas::Field horSpecVarField("horSpecVar", make_datatype<double>(),
        make_shape(trans_->nw(), nz, 1));

      // Get horizontal spectral variance view
      auto horSpecVarView = make_view<double, 3>(horSpecVarField);

      // Set horizontal spectral variance to zero
      horSpecVarView.assign(0.0);

      // Compute horizontal spectral variance
      for (size_t js = 0; js < trans_->ns(); ++js) {
        if (trans_->jq(js) == 0) {
          const size_t jw = trans_->jw(js);
          for (size_t jzI = 0; jzI < nz; ++jzI) {
            horSpecVarView(jw, jzI, 0) += horCorSpView(js, jzI);
          }
        }
      }

      // Reduce horizontal spectral variance
      trans_->reduceCov(horSpecVarField);

      // Compute horizontal spectral standard-deviation
      for (size_t jw = 0; jw < trans_->nw(); ++jw) {
        for (size_t jzI = 0; jzI < nz; ++jzI) {
          horSpecVarView(jw, jzI, 0) *= trans_->spNormSumInv(jw);
          horSpecVarView(jw, jzI, 0) = std::max(horSpecVarView(jw, jzI, 0), 0.0);
          horSpecVarView(jw, jzI, 0) = std::sqrt(horSpecVarView(jw, jzI, 0));
        }
      }

      // Compute correlation square-root
      for (size_t jw = 0; jw < trans_->nw(); ++jw) {
        for (size_t jzI = 0; jzI < nz; ++jzI) {
          for (size_t jzJ = 0; jzJ < nz; ++jzJ) {
            corSqrtView(jw, jzI, jzJ) = horSpecVarView(jw, jzI, 0)*es.eigenvectors().col(jzJ)[jzI]
              *std::sqrt(es.eigenvalues()[jzJ]);
          }
        }
      }

      // Create covariance field
      createField3D("cov", trans_->nw(), var, data_);

      // Get covariance view
      auto covView = getView3D("cov", var, data_);

      // Compute covariance matrix
      for (size_t jw = 0; jw < trans_->nw(); ++jw) {
        for (size_t jzI = 0; jzI < nz; ++jzI) {
          for (size_t jzJ = 0; jzJ < nz; ++jzJ) {
            for (size_t jz3 = 0; jz3 < nz; ++jz3) {
              covView(jw, jzI, jzJ) += corSqrtView(jw, jzI, jz3)*corSqrtView(jw, jzJ, jz3);
            }
          }
        }
      }
    }
  } else {
    // Ensemble-based calibration
    ASSERT(ne > 2);

    for (const auto & var : activeVars_) {
      // Get number of levels
      const size_t nz = var.getLevels();

      // Create covariance field
      createField3D("cov", trans_->nw(), var, data_);

      // Get covariance view
      auto covView = getView3D("cov", var, data_);

      // Loop over ensemble members
      for (size_t je = 0; je < ne; ++je) {
        // Get view
        const auto view = getView2D(var, fsetEns[je]);

        // Update covariance (lower triangle)
        for (size_t js = 0; js < trans_->ns(); ++js) {
          for (size_t jw = 0; jw < trans_->nw(); ++jw) {
            if (trans_->includeWavenumber(js, jw)) {
              const double factor = trans_->spNorm(js);
              for (size_t jzI = 0; jzI < nz; ++jzI) {
                for (size_t jzJ = 0; jzJ < jzI+1; ++jzJ) {
                  covView(jw, jzI, jzJ) += factor*view(js, jzJ)*view(js, jzI);
                }
              }
            }
          }
        }
      }

      // Transpose lower triangle
      for (size_t jw = 0; jw < trans_->nw(); ++jw) {
        for (size_t jzI = 0; jzI < nz; ++jzI) {
          for (size_t jzJ = 0; jzJ < jzI; ++jzJ) {
            covView(jw, jzJ, jzI) = covView(jw, jzI, jzJ);
          }
        }
      }

      // Get covariance field
      auto covField = getField("cov", var, data_);

      // Reduce and normalize covariance
      trans_->reduceNormalizeCov(ne-1, covField);

      // Filter covariance
      trans_->filterCov(Lf_, covField);
    }
  }

  // Compute square-root
  computeSquareRoot();

  if (ne == 0) {
    // Update the standard-deviation
    for (const auto & var : activeVars_) {
      // Get number of levels
      const size_t nz = var.getLevels();

      // Get standard-deviation profile if present (default = 1.0)
      std::vector<double> stdDev(nz, 1.0);
      for (const auto & profile : params_.calibration.value()->profiles.value()) {
        if (profile.variable.value() == var.name()) {
          if (profile.stdDev.value() != boost::none) {
            stdDev = *profile.stdDev.value();
          }
        }
      }

      // Get standard-deviation view
      auto stdDevView = getViewProfile("stdDev", var, data_);

      // Copy standard-deviation
      for (size_t jzI = 0; jzI < nz; ++jzI) {
        stdDevView(jzI) = stdDev[jzI];
      }
    }
  }

  // Print norms
  print(oops::Log::test());

  oops::Log::trace() << classname() << "::directCalibration done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierCovariance::iterativeCalibrationInit() {
  oops::Log::trace() << classname() << "::iterativeCalibrationInit starting" << std::endl;

  // Initialize iterative counters with zeroes
  iterativeN_ = 0;

  for (const auto & var : activeVars_) {
    // Create perturbation field
    createField2D("pert", trans_->ns(), var, data_);

    // Create mean field
    createField2D("mean", trans_->ns(), var, data_);

    // Create covariance field
    createField3D("cov", trans_->nw(), var, data_);
  }

  oops::Log::trace() << classname() << "::iterativeCalibrationInit done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierCovariance::iterativeCalibrationUpdate(const oops::FieldSet3D & fset) {
  oops::Log::trace() << classname() << "::iterativeCalibrationUpdate starting" << std::endl;

  // Increment ensemble index (ie = ie + 1)
  ++iterativeN_;

  // Sub-ensemble index
  const size_t ie = (params_.calibration.value()->subEnsSize.value() > 0) ?
    ((iterativeN_-1)%params_.calibration.value()->subEnsSize.value())+1 : iterativeN_;

  for (const auto & var : activeVars_) {
    // Get number of output levels
    const size_t nz = var.getLevels();

    // Get member view
    const auto view = getView2D(var, fset);

    // Get perturbation view
    auto pertView = getView2D("pert", var, data_);

    // Get mean view
    auto meanView = getView2D("mean", var, data_);

    if (ie == 1) {
      // Reset mean if a new sub-ensemble starts
      meanView.assign(0.0);
    }

    // Remove mean (pert = state - mean)
    for (size_t js = 0; js < trans_->ns(); ++js) {
      for (size_t jz = 0; jz < nz; ++jz) {
        pertView(js, jz) = view(js, jz) - meanView(js, jz);
      }
    }

    if (ie > 1) {
      // Get covariance view
      auto covView = getView3D("cov", var, data_);

      // Update covariance (cov = cov + (ie-1)/ie * pert1 * pert2)
      for (size_t js = 0; js < trans_->ns(); ++js) {
        for (size_t jw = 0; jw < trans_->nw(); ++jw) {
          if (trans_->includeWavenumber(js, jw)) {
            const double factor = trans_->spNorm(js)
              *static_cast<double>(ie-1)/static_cast<double>(ie);
            for (size_t jzI = 0; jzI < nz; ++jzI) {
              for (size_t jzJ = 0; jzJ < nz; ++jzJ) {
                covView(jw, jzI, jzJ) += factor*pertView(js, jzJ)*pertView(js, jzI);
              }
            }
          }
        }
      }
    }

    // Update mean (mean = mean + 1 / ie * pert)
    const double factor = 1.0/static_cast<double>(ie);
    for (size_t js = 0; js < trans_->ns(); ++js) {
      for (size_t jz = 0; jz < nz; ++jz) {
        meanView(js, jz) += factor*pertView(js, jz);
      }
    }
  }

  oops::Log::trace() << classname() << "::iterativeCalibrationUpdate done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierCovariance::iterativeCalibrationFinal() {
  oops::Log::trace() << classname() << "::iterativeCalibrationFinal starting" << std::endl;

  // Compute number of sub-ensembles
  size_t nSubEns = 1;
  if (params_.calibration.value()->subEnsSize.value() > 0) {
    // Check sub-ensembles size
    ASSERT(iterativeN_%params_.calibration.value()->subEnsSize.value() == 0);

    // Get number of sub-ensembles
    nSubEns = iterativeN_/params_.calibration.value()->subEnsSize.value();
  }

  for (const auto & var : activeVars_) {
    // Get covariance field
    auto covField = getField("cov", var, data_);

    // Reduce and normalize covariance
    trans_->reduceNormalizeCov(iterativeN_-nSubEns, covField);

    // Filter covariance
    trans_->filterCov(Lf_, covField);
  }

  // Compute square-root
  computeSquareRoot();

  // Print norms
  print(oops::Log::test());

  oops::Log::trace() << classname() << "::iterativeCalibrationFinal done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierCovariance::write() const {
  oops::Log::trace() << classname() << "::write starting" << std::endl;

  if (params_.write.value() != boost::none) {
    // Create covariance fieldset
    atlas::FieldSet covData;

    if (params_.write.value()->writeCovariance.value()) {
      // Compute covariance from correlation square-root and standard-deviation if it is missing
      computeCovariance(covData);
    }

    // NetCDF IDs
    int retval, ncId, nwId, nzIId, nzJId, dCorSqrtId[3], dStdDevId[1], dCovId[3],
      corSqrtId[activeVars_.size()], stdDevId[activeVars_.size()], covId[activeVars_.size()];

    // NetCDF file path
    const std::string ncFilePath = params_.write.value()->outputFile.value();

    // Definition mode
    size_t jvar = 0;

    if (comm_.rank() == 0) {
      // Create NetCDF file
      if ((retval = nc_create(ncFilePath.c_str(), NC_64BIT_OFFSET | NC_CLOBBER, &ncId)))
        ERR(retval, ncFilePath);

      // Create horizontal dimension
      if ((retval = nc_def_dim(ncId, "nwGlb", trans_->nwGlb(), &nwId)))
        ERR(retval, "nwGlb");

      // Dimensions arrays, horizontal part
      dCorSqrtId[0] = nwId;
      dCovId[0] = nwId;

      for (const auto & var : activeVars_) {
        // Get number of levels
        const size_t nz = var.getLevels();

        // Create vertical dimensions
        const std::string nzIName = "nzI_" + var.name();
        const std::string nzJName = "nzJ_" + var.name();
        if ((retval = nc_def_dim(ncId, nzIName.c_str(), nz, &nzIId))) ERR(retval, nzIName);
        if ((retval = nc_def_dim(ncId, nzJName.c_str(), nz, &nzJId))) ERR(retval, nzJName);

        // Dimensions array, vertical part
        dCorSqrtId[1] = nzIId;
        dCorSqrtId[2] = nzJId;
        dStdDevId[0] = nzIId;
        dCovId[1] = nzIId;
        dCovId[2] = nzJId;

        // Get correlation square-root field name
        const std::string corSqrtFieldName = fieldName("corSqrt", var);

        // Get standard-deviation field name
        const std::string stdDevFieldName = fieldName("stdDev", var);

        // Get correlation square-root field name
        const std::string covFieldName = fieldName("cov", var);

        // Define variables
        if ((retval = nc_def_var(ncId, corSqrtFieldName.c_str(), NC_DOUBLE, 3, dCorSqrtId,
          &corSqrtId[jvar]))) ERR(retval, corSqrtFieldName);
        if ((retval = nc_def_var(ncId, stdDevFieldName.c_str(), NC_DOUBLE, 1, dStdDevId,
          &stdDevId[jvar]))) ERR(retval, stdDevFieldName);
        if (params_.write.value()->writeCovariance.value()) {
          if ((retval = nc_def_var(ncId, covFieldName.c_str(), NC_DOUBLE, 3, dCovId,
            &covId[jvar]))) ERR(retval, covFieldName);
        }

        // Update index
        ++jvar;
      }

      // End definition mode
      if ((retval = nc_enddef(ncId))) ERR(retval, ncFilePath);
    }

    // Data mode
    jvar = 0;

    for (const auto & var : activeVars_) {
      // Get number of levels
      const size_t nz = var.getLevels();

      // Define global vectors
      std::vector<double> corSqrtVecGlb;
      std::vector<double> stdDevVec;
      std::vector<double> covVecGlb;

      // Get correlation square-root field
      const auto corSqrtField = getField("corSqrt", var, data_);

      // Gather correlation square-root vector
      trans_->gatherCov(corSqrtField, corSqrtVecGlb);

      if (params_.write.value()->writeCovariance.value()) {
        // Get correlation square-root field
        const auto covField = getField("cov", var, covData);

        // Gather correlation square-root vector
        trans_->gatherCov(covField, covVecGlb);
      }

      if (comm_.rank() == 0) {
        // Get standard-deviation view
        const auto stdDevView = getViewProfile("stdDev", var, data_);

        // Allocate global standard-deviation vector
        stdDevVec.resize(nz);

        // Serialize
        for (size_t jz = 0; jz < nz; ++jz) {
          stdDevVec[jz] = stdDevView(jz);
        }

        // Write data
        if ((retval = nc_put_var_double(ncId, corSqrtId[jvar], corSqrtVecGlb.data())))
          ERR(retval, var.name() + "_corSqrt");
        if ((retval = nc_put_var_double(ncId, stdDevId[jvar], stdDevVec.data())))
          ERR(retval, var.name() + "_stdDev");
        if (params_.write.value()->writeCovariance.value()) {
          if ((retval = nc_put_var_double(ncId, covId[jvar], covVecGlb.data())))
            ERR(retval, var.name() + "_cov");
        }
        ++jvar;
      }
    }

    if (comm_.rank() == 0) {
      // Close file
      if ((retval = nc_close(ncId))) ERR(retval, ncFilePath);
    }
  }
}

// -----------------------------------------------------------------------------

void BifourierCovariance::readCovariance() {
  oops::Log::trace() << classname() << "::readCovariance starting" << std::endl;

  for (const auto & var : activeVars_) {
    // Create covariance field
    createField3D("oldCov", trans_->nw(), var, data_);
  }

  // NetCDF file path
  const std::string ncFilePath = *params_.calibration.value()->oldCovInputFile.value();

  // NetCDF IDs
  int ncId, retval, nwId, nzIId, nzJId, covId;
  size_t nwGlbFromFile, nzIFromFile, nzJFromFile;

  if (comm_.rank() == 0) {
    // Open NetCDF file
    if ((retval = nc_open(ncFilePath.c_str(), NC_NOWRITE, &ncId))) ERR(retval, ncFilePath);
  }

  for (const auto & var : activeVars_) {
    // Get number of levels
    const size_t nz = var.getLevels();

    // Define global covariance vector
    std::vector<double> covVecGlb;

    if (comm_.rank() == 0) {
      // Check dimensions
      const std::string nzIName = "nzI_" + var.name();
      const std::string nzJName = "nzJ_" + var.name();
      if ((retval = nc_inq_dimid(ncId, "nwGlb", &nwId))) ERR(retval, "nwGlb");
      if ((retval = nc_inq_dimid(ncId, nzIName.c_str(), &nzIId))) ERR(retval, nzIName);
      if ((retval = nc_inq_dimid(ncId, nzJName.c_str(), &nzJId))) ERR(retval, nzJName);
      if ((retval = nc_inq_dimlen(ncId, nwId, &nwGlbFromFile))) ERR(retval, "nwGlb");
      if ((retval = nc_inq_dimlen(ncId, nzIId, &nzIFromFile))) ERR(retval, nzIName);
      if ((retval = nc_inq_dimlen(ncId, nzJId, &nzJFromFile))) ERR(retval, nzJName);
      ASSERT(nwGlbFromFile == trans_->nwGlb());
      ASSERT(nzIFromFile == nz);
      ASSERT(nzJFromFile == nz);

      // Get covariance field name
      const std::string covFieldName = fieldName("cov", var);

      // Get variables ID
      if ((retval = nc_inq_varid(ncId, covFieldName.c_str(), &covId)))
        ERR(retval, covFieldName);

      // Allocate global covariance vector
      covVecGlb.resize(trans_->nwGlb()*nz*nz);

      // Read data
      if ((retval = nc_get_var_double(ncId, covId, covVecGlb.data())))
        ERR(retval, covFieldName);
    }

    // Get covariance field
    auto covField = getField("oldCov", var, data_);

    // Scatter covariance vector
    trans_->scatterCov(covVecGlb, covField);
  }

  if (comm_.rank() == 0) {
    // Close file
    if ((retval = nc_close(ncId))) ERR(retval, ncFilePath);
  }

  oops::Log::trace() << classname() << "::readCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierCovariance::computeSquareRoot() {
  oops::Log::trace() << classname() << "::computeSquareRoot starting" << std::endl;

  // Combine covariance with an old covariance
  if (params_.calibration.value() != boost::none) {
    if (params_.calibration.value()->oldCovInputFile.value() != boost::none) {
      // Check parameters consistency
      ASSERT(params_.calibration.value()->halfLife.value() != boost::none);

      // Read old covariance
      readCovariance();

      // Define update factor
      const double halfLife = *params_.calibration.value()->halfLife.value();
      const double alphaInf = 1.0-std::exp(-std::log(2.0)/halfLife);
      double updateFactor = alphaInf;
      if (params_.calibration.value()->cycleIndex.value() != boost::none) {
        const size_t cycleIndex = *params_.calibration.value()->cycleIndex.value();
        ASSERT(cycleIndex > 0);
        updateFactor /= 1.0-std::pow(1.0-alphaInf, static_cast<double>(cycleIndex+1));
      }

      for (const auto & var : activeVars_) {
        // Get number of levels
        const size_t nz = var.getLevels();

        // Get old covariance view
        const auto oldCovView = getView3D("oldCov", var, data_);

        // Get covariance view
        auto covView = getView3D("cov", var, data_);

        // Combine covariances
        for (size_t jw = 0; jw < trans_->nw(); ++jw) {
          for (size_t jzI = 0; jzI < nz; ++jzI) {
            for (size_t jzJ = 0; jzJ < nz; ++jzJ) {
              covView(jw, jzI, jzJ) = updateFactor*covView(jw, jzI, jzJ)
                + (1.0-updateFactor)*oldCovView(jw, jzI, jzJ);
            }
          }
        }
      }
    }
  }

  for (const auto & var : activeVars_) {
    // Get number of levels
    const size_t nz = var.getLevels();

    // Create correlation square-root field
    createField3D("corSqrt", trans_->nw(), var, data_);

    // Create standard-deviation field
    createFieldProfile("stdDev", var, data_);

    // Get covariance view
    const auto covView = getView3D("cov", var, data_);

    // Get correlation square-root view
    auto corSqrtView = getView3D("corSqrt", var, data_);

    // Get standard-deviation view
    auto stdDevView = getViewProfile("stdDev", var, data_);

    // Compute vertical correlation square-root
    const double zeps = 1.0e-99;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    for (size_t jw = 0; jw < trans_->nw(); ++jw) {
      // Compute vertical correlation matrix
      Eigen::MatrixXd vertCor(nz, nz);
      for (size_t jz2 = 0; jz2 < nz; ++jz2) {
        for (size_t jz1 = 0; jz1 < nz; ++jz1) {
          vertCor(jz2, jz1) = covView(jw, jz2, jz1) /
            std::sqrt(std::max(zeps, covView(jw, jz1, jz1))*std::max(zeps, covView(jw, jz2, jz2)));
        }
      }

      // Compute eigendecomposition
      es.compute(vertCor);

      // Store covariance square-root
      for (size_t jz2 = 0; jz2 < nz; ++jz2) {
        for (size_t jz1 = 0; jz1 < nz; ++jz1) {
          corSqrtView(jw, jz2, jz1) = es.eigenvectors().col(jz1)[jz2]
            *std::sqrt(es.eigenvalues()[jz1]);
        }
      }
    }

    // Compute sum over wavenumbers
    std::vector<double> sum(nz, 0.0);
    for (size_t jw = 0; jw < trans_->nwRoot(); ++jw) {
      for (size_t jz = 0; jz < nz; ++jz) {
        sum[jz] += covView(jw, jz, jz);
      }
    }

    // Communication
    comm_.allReduceInPlace(sum.begin(), sum.end(), eckit::mpi::sum());

    // Create horizontal covariance field
    atlas::Field horCovField("horCov", make_datatype<double>(),
      make_shape(trans_->nw(), nz));

    // Get horizontal covariance view
    auto horCovView = make_view<double, 2>(horCovField);

    // Normalize with the sum
    for (size_t jw = 0; jw < trans_->nw(); ++jw) {
      for (size_t jz = 0; jz < nz; ++jz) {
        horCovView(jw, jz) = covView(jw, jz, jz)/sum[jz];
      }
    }

    // Take square-root
    for (size_t jw = 0; jw < trans_->nw(); ++jw) {
      for (size_t jz = 0; jz < nz; ++jz) {
        horCovView(jw, jz) = std::sqrt(horCovView(jw, jz));
      }
    }

    // Create standard-deviation profiles
    std::vector<double> vertStd(nz, 0.0);

    // Compute vertical variance for each level
    for (size_t jw = 0; jw < trans_->nwRoot(); ++jw) {
      for (size_t jz = 0; jz < nz; ++jz) {
        vertStd[jz] += covView(jw, jz, jz);
      }
    }

    // Communication
    comm_.allReduceInPlace(vertStd.begin(), vertStd.end(), eckit::mpi::sum());

    // Take variance square-root
    for (size_t jz = 0; jz < nz; ++jz) {
      vertStd[jz] = std::sqrt(vertStd[jz]);
    }

    // Merge contributions
    for (size_t jw = 0; jw < trans_->nw(); ++jw) {
      for (size_t jz2 = 0; jz2 < nz; ++jz2) {
        for (size_t jz1 = 0; jz1 < nz; ++jz1) {
          corSqrtView(jw, jz2, jz1) *= horCovView(jw, jz2)*vertStd[jz2];
        }
      }
    }

    // Compute variance
    std::vector<double> variance(nz, 0.0);
    for (size_t jz2 = 0; jz2 < nz; ++jz2) {
      for (size_t js = 0; js < trans_->ns(); ++js) {
        if (trans_->jq(js) == 0) {
          const size_t jw = trans_->jw(js);
          for (size_t jz1 = 0; jz1 < nz; ++jz1) {
            variance[jz2] += corSqrtView(jw, jz2, jz1)*corSqrtView(jw, jz2, jz1)
              *trans_->spNorm(js)*trans_->spNorm(js);
          }
        }
      }
    }

    // Communication
    comm_.allReduceInPlace(variance.begin(), variance.end(), eckit::mpi::sum());

    // Compute standard-deviation
    for (size_t jz = 0; jz < nz; ++jz) {
      stdDevView(jz) = std::sqrt(variance[jz]);
    }

    // Apply inverse standard-deviation to normalize correlation square-root
    for (size_t jz2 = 0; jz2 < nz; ++jz2) {
      ASSERT(stdDevView(jz2) > 0.0);
      const double norm = 1.0/stdDevView(jz2);
      for (size_t jw = 0; jw < trans_->nw(); ++jw) {
        for (size_t jz1 = 0; jz1 < nz; ++jz1) {
          corSqrtView(jw, jz2, jz1) *= norm;
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::computeSquareRoot done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierCovariance::computeCovariance(atlas::FieldSet & covData) const {
  oops::Log::trace() << classname() << "::computeCovariance starting" << std::endl;

  for (const auto & var : activeVars_) {
    // Get covariance field name
    const auto covFieldName = fieldName("cov", var);

    if (data_.has(covFieldName)) {
      // Share pointer
      covData.add(data_[covFieldName]);
    } else {
      // Get number of levels
      const size_t nz = var.getLevels();

      // Create covariance field
      createField3D("cov", trans_->nw(), var, covData);

      // Get covariance view
      auto covView = getView3D("cov", var, covData);

      // Get correlation square-root view
      const auto corSqrtView = getView3D("corSqrt", var, data_);

      // Get standard-deviation view
      const auto stdDevView = getViewProfile("stdDev", var, data_);

      // Compute covariance matrix
      for (size_t jw = 0; jw < trans_->nw(); ++jw) {
        for (size_t jzI = 0; jzI < nz; ++jzI) {
          for (size_t jzJ = 0; jzJ < nz; ++jzJ) {
            // Compute correlation
            for (size_t jz3 = 0; jz3 < nz; ++jz3) {
              covView(jw, jzI, jzJ) += corSqrtView(jw, jzI, jz3)*corSqrtView(jw, jzJ, jz3);
            }

            // Multiply by standard-deviations
            covView(jw, jzI, jzJ) *= stdDevView(jzJ)*stdDevView(jzI);
          }
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::computeCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierCovariance::print(std::ostream & os) const {
  // Print norms
  os << "Covariance norms: " << std::endl;
  for (const auto & var : activeVars_) {
    // Get number of levels
    const size_t nz = var.getLevels();

    // Get correlation square-root field
    const auto corSqrtField = getField("corSqrt", var, data_);

    // Get standard-deviation view
    const auto stdDevView = getViewProfile("stdDev", var, data_);

    // Compute standard deviation squared norm
    double stdDevNorm = 0.0;
    for (size_t jz = 0; jz < nz; ++jz) {
      stdDevNorm += stdDevView(jz)*stdDevView(jz);
    }

    // Compute standard deviation norm
    stdDevNorm = std::sqrt(stdDevNorm);

    // Print norms
    os << "- " << var.name() << ": " << std::endl;
    os << "  + correlation square-root: " << trans_->normCov(corSqrtField) << std::endl;
    os << "  + standard-deviation:      " << stdDevNorm << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber

