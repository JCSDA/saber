/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/BifourierBalance.h"

#include <Eigen/Dense>
#include <netcdf.h>

#include <algorithm>

#include "oops/util/FieldSetOperations.h"

#include "saber/bifourier/BifourierUtilities.h"

#define ERR(e, msg) {std::string s(nc_strerror(e)); \
  throw eckit::Exception(s + " : " + msg, Here());}

using atlas::array::make_datatype;
using atlas::array::make_shape;
using atlas::array::make_view;

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

static OuterBlockMaker<BifourierBalance> makerBifourierBalance_("BifourierBalance");

// -----------------------------------------------------------------------------

BifourierBalance::BifourierBalance(const oops::GeometryData & outerGeometryData,
                                   const oops::Variables & outerVars,
                                   const eckit::Configuration & covarConfig,
                                   const Parameters_ & params,
                                   const oops::FieldSet3D & xb,
                                   const oops::FieldSet3D & fg)
  : OuterBlockBase(params, xb.validTime(), outerGeometryData, outerVars),
    innerGeometryData_(outerGeometryData),
    comm_(outerGeometryData.comm()),
    innerVars_(outerVars),
    params_(params),
    Lf_(params_.calibration.value() ? params_.calibration.value()->filteringScale.value() : 0),
    trans_(transStore_.retrieveTransform(outerGeometryData))
{
  oops::Log::trace() << classname() << "::BifourierBalance starting" << std::endl;

  // Initialize components counter
  nCmp_ = 0;

  // Check the rows consistency: blocks should be balance to create a lower triangular matrix
  for (const auto & row : params_.rows.value()) {
    // Get output variable name
    const std::string outputVarName = row.outputVar.value();

    // Check that output variable is present in block variables
    ASSERT(innerVars_.has(outputVarName));

    for (const auto & inputVarName : row.inputVars.value()) {
      // Check that input variable is present in block variables
      ASSERT(innerVars_.has(inputVarName));

      // Input variable should be already present in the balance variables
      ASSERT(balVars_.has(inputVarName));

      // Update counter
      ++nCmp_;
    }

    // Add output variable into balance variables
    balVars_.push_back(innerVars_[outputVarName]);
  }

  // Check that all variables are present (but maybe in a different order)
  ASSERT(balVars_ <= innerVars_);

  // Extra variable for auto-covariance + balance variables
  if (params_.extraVar.value()) {
    ASSERT(!balVars_.has(*params_.extraVar.value()));
    balVarsExt_.push_back(innerVars_[*params_.extraVar.value()]);
  }
  for (const auto & balVar : balVars_) {
    balVarsExt_.push_back(balVar);
  }

  oops::Log::trace() << classname() << "::BifourierBalance done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierBalance::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  for (int jr = params_.rows.value().size()-1; jr >= 0; --jr) {
    // Get row
    const auto & row = params_.rows.value()[jr];

    // Get output variable
    const oops::Variable outputVar = balVars_[row.outputVar.value()];

    // Get number of output levels
    const size_t nzI = outputVar.getLevels();

    // Get output view
    auto outputView = getView2D(outputVar, fset);

    for (const auto & inputVarName : row.inputVars.value()) {
      // Get input variable
      const oops::Variable inputVar = balVars_[inputVarName];

      // Get number of input levels
      const size_t nzJ = inputVar.getLevels();

      // Get number of input levels and input view
      const auto inputView = getView2D(inputVar, fset);

      // Get regression view
      const auto regView = getView3D("reg", outputVar, inputVar, data_);

      // Apply matrix
      for (size_t js = 0; js < trans_->ns(); ++js) {
        const size_t jw = trans_->jw(js);
        for (size_t jzI = 0; jzI < nzI; ++jzI) {
          for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
            outputView(js, jzI) += regView(jw, jzI, jzJ)*inputView(js, jzJ);
          }
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierBalance::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  for (const auto & row : params_.rows.value()) {
    // Get output variable
    const oops::Variable outputVar = balVars_[row.outputVar.value()];

    // Get number of output levels
    const size_t nzI = outputVar.getLevels();

    // Get output view
    const auto outputView = getView2D(outputVar, fset);

    for (const auto & inputVarName : row.inputVars.value()) {
      // Get input variable
      const oops::Variable inputVar = balVars_[inputVarName];

      // Get number of input levels
      const size_t nzJ = inputVar.getLevels();

      // Get number of input levels and input view
      auto inputView = getView2D(inputVar, fset);

      // Get regression view
      const auto regView = getView3D("reg", outputVar, inputVar, data_);

      // Apply matrix
      for (size_t js = 0; js < trans_->ns(); ++js) {
        const size_t jw = trans_->jw(js);
        for (size_t jzI = 0; jzI < nzI; ++jzI) {
          for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
            inputView(js, jzJ) += regView(jw, jzI, jzJ)*outputView(js, jzI);
          }
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierBalance::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;

  for (const auto & row : params_.rows.value()) {
    // Get output variable
    const oops::Variable outputVar = balVars_[row.outputVar.value()];

    // Get number of output levels
    const size_t nzI = outputVar.getLevels();

    // Get output view
    auto outputView = getView2D(outputVar, fset);;

    for (const auto & inputVarName : row.inputVars.value()) {
      // Get input variable
      const oops::Variable inputVar = balVars_[inputVarName];

      // Get number of input levels
      const size_t nzJ = inputVar.getLevels();

      // Get number of input levels and input view
      const auto inputView = getView2D(inputVar, fset);

      // Get regression view
      const auto regView = getView3D("reg", outputVar, inputVar, data_);

      // Apply matrix
      for (size_t js = 0; js < trans_->ns(); ++js) {
        const size_t jw = trans_->jw(js);
        for (size_t jzI = 0; jzI < nzI; ++jzI) {
          for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
            outputView(js, jzI) -= regView(jw, jzI, jzJ)*inputView(js, jzJ);
          }
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierBalance::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;

  for (const auto & row : params_.rows.value()) {
    // Get output variable
    const oops::Variable outputVar = balVars_[row.outputVar.value()];
    for (const auto & inputVarName : row.inputVars.value()) {
      // Get input variable
      const oops::Variable inputVar = balVars_[inputVarName];

      // Create regression field
      createField3D("reg", trans_->nw(), outputVar, inputVar, data_);
    }
  }

  // NetCDF file path
  const std::string ncFilePath = params_.read.value()->inputFile.value();

  // NetCDF IDs
  int ncId, retval, nwGlbId, nzIId, nzJId, varId;
  size_t nwGlbFromFile, nzIFromFile, nzJFromFile;

  if (comm_.rank() == 0) {
    // Open NetCDF file
    if ((retval = nc_open(ncFilePath.c_str(), NC_NOWRITE, &ncId))) ERR(retval, ncFilePath);
  }

  for (const auto & row : params_.rows.value()) {
    // Get output variable
    const oops::Variable outputVar = balVars_[row.outputVar.value()];

    // Get number of output levels
    const size_t nzI = outputVar.getLevels();

    // Define global regression vector
    std::vector<double> regVecGlb;

    for (const auto & inputVarName : row.inputVars.value()) {
      // Get input variable
      const oops::Variable inputVar = balVars_[inputVarName];

      // Get number of input levels
      const size_t nzJ = inputVar.getLevels();

      if (comm_.rank() == 0) {
        // Check dimensions
        const std::string nzIName = "nzI_" + outputVar.name();
        const std::string nzJName = "nzJ_" + inputVar.name();
        if ((retval = nc_inq_dimid(ncId, "nwGlb", &nwGlbId))) ERR(retval, "nwGlb");
        if ((retval = nc_inq_dimid(ncId, nzIName.c_str(), &nzIId))) ERR(retval, nzIName);
        if ((retval = nc_inq_dimid(ncId, nzJName.c_str(), &nzJId))) ERR(retval, nzJName);
        if ((retval = nc_inq_dimlen(ncId, nwGlbId, &nwGlbFromFile))) ERR(retval, "nwGlb");
        if ((retval = nc_inq_dimlen(ncId, nzIId, &nzIFromFile))) ERR(retval, nzIName);
        if ((retval = nc_inq_dimlen(ncId, nzJId, &nzJFromFile))) ERR(retval, nzJName);
        ASSERT(nwGlbFromFile == trans_->nwGlb());
        ASSERT(nzIFromFile == nzI);
        ASSERT(nzJFromFile == nzJ);

        // Allocate global regression vector
        regVecGlb.resize(trans_->nwGlb()*nzI*nzJ);

        // Read data
        const std::string regFieldName = fieldName("reg", outputVar, inputVar);
        if ((retval = nc_inq_varid(ncId, regFieldName.c_str(), &varId)))
          ERR(retval, regFieldName);
        if ((retval = nc_get_var_double(ncId, varId, regVecGlb.data()))) ERR(retval, regFieldName);
      }

      // Get regression field
      auto regField = getField("reg", outputVar, inputVar, data_);

      // Scatter regression vector
      trans_->scatterCov(regVecGlb, regField);
    }
  }

  if (comm_.rank() == 0) {
    // Close file
    if ((retval = nc_close(ncId))) ERR(retval, ncFilePath);
  }

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierBalance::directCalibration(const oops::FieldSets & fsetEns) {
  oops::Log::trace() << classname() << "::directCalibration starting" << std::endl;

  // Check ensemble size
  const size_t ne = fsetEns.ens_size();
  ASSERT(ne > 2);

  // Check if xx-covariances must be computed
  bool computeXxCov = params_.calibration.value()->fullRecursiveInverse.value();
  if (params_.write.value()) {
    computeXxCov = computeXxCov || params_.write.value()->writeCovariance.value();
  }

  if (computeXxCov) {
    // Estimate xx-covariance
    for (const auto & outputVar : balVarsExt_) {
      // Get number of output levels
      const size_t nzI = outputVar.getLevels();

      for (const auto & inputVar : xxCovVars(outputVar)) {
        // Get number of input levels
        const size_t nzJ = inputVar.getLevels();

        // Create xx-covariance field
        createField3D("xxCov", trans_->nw(), outputVar, inputVar, data_);

        // Get xx-covariance view
        auto xxCovView = getView3D("xxCov", outputVar, inputVar, data_);

        // Loop over ensemble members
        for (size_t je = 0; je < ne; ++je) {
          // Get member views
          const auto outputView = getView2D(outputVar, fsetEns[je]);
          const auto inputView = getView2D(inputVar, fsetEns[je]);

          // Update xx-covariance
          for (size_t js = 0; js < trans_->ns(); ++js) {
            for (size_t jw = 0; jw < trans_->nw(); ++jw) {
              if (trans_->includeWavenumber(js, jw)) {
                const double factor = trans_->spNorm(js);
                for (size_t jzI = 0; jzI < nzI; ++jzI) {
                  for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
                    xxCovView(jw, jzI, jzJ) += factor*inputView(js, jzJ)*outputView(js, jzI);
                  }
                }
              }
            }
          }
        }

        // Get xx-covariance field
        auto xxCovField = getField("xxCov", outputVar, inputVar, data_);

        // Reduce and normalize xx-covariance
        trans_->reduceNormalizeCov(ne-1, xxCovField);
      }
    }
  }

  if (params_.calibration.value()->fullRecursiveInverse.value()) {
    // Using the full recursive inverse formula

    // Compute regressions from xx-covariances
    computeRegressionsFromCovariances();
  } else {
    // Using the partial recursive inverse formula

    // Copy ensemble
    oops::FieldSets fsetEnsUnbal(fsetEns);

    for (const auto & row : params_.rows.value()) {
      // Get output variable
      const oops::Variable outputVar = balVars_[row.outputVar.value()];

      // Get number of output levels
      const size_t nzI = outputVar.getLevels();

      for (const auto & inputVarName : row.inputVars.value()) {
        // Get input variable
        const oops::Variable inputVar = balVars_[inputVarName];

        // Get number of input levels
        const size_t nzJ = inputVar.getLevels();

        // Create xv-covariance field
        createField3D("xvCov", trans_->nw(), outputVar, inputVar, data_);

        // Get xv-covariance view
        auto xvCovView = getView3D("xvCov", outputVar, inputVar, data_);

        // Loop over ensemble members
        for (size_t je = 0; je < ne; ++je) {
          // Get member views
          const auto outputView = getView2D(outputVar, fsetEns[je]);
          const auto inputUnbalView = getView2D(inputVar, fsetEnsUnbal[je]);

          // Update xv-covariance
          for (size_t js = 0; js < trans_->ns(); ++js) {
            for (size_t jw = 0; jw < trans_->nw(); ++jw) {
              if (trans_->includeWavenumber(js, jw)) {
                const double factor = trans_->spNorm(js);
                for (size_t jzI = 0; jzI < nzI; ++jzI) {
                  for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
                    xvCovView(jw, jzI, jzJ) += factor*inputUnbalView(js, jzJ)
                      *outputView(js, jzI);
                  }
                }
              }
            }
          }
        }

        // Get xv-covariance field
        auto xvCovField = getField("xvCov", outputVar, inputVar, data_);

        // Reduce and normalize xv-covariance
        trans_->reduceNormalizeCov(ne-1, xvCovField);

        // Filter xv-covariance
        trans_->filterCov(Lf_, xvCovField);
      }

      // Compute regression
      computeRegression(row.inputVars.value(), outputVar);

      // Update ensemble members
      for (const auto & inputVarName : row.inputVars.value()) {
        // Get input variable
        const oops::Variable inputVar = balVars_[inputVarName];

        // Get number of input levels
        const size_t nzJ = inputVar.getLevels();

        // Get regression view
        const auto regView = getView3D("reg", outputVar, inputVar, data_);

        for (size_t je = 0; je < ne; ++je) {
          // Get member views
          const auto inputUnbalView = getView2D(inputVar, fsetEnsUnbal[je]);
          auto outputUnbalView = getView2D(outputVar, fsetEnsUnbal[je]);

          // Unbalance perturbation
          for (size_t js = 0; js < trans_->ns(); ++js) {
            const size_t jw = trans_->jw(js);
            for (size_t jzI = 0; jzI < nzI; ++jzI) {
              for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
                outputUnbalView(js, jzI) -= regView(jw, jzI, jzJ)*inputUnbalView(js, jzJ);
              }
            }
          }
        }
      }

      // vv-covariance variables
      std::vector<std::string> vvCovVars = row.inputVars.value();
      vvCovVars.push_back(outputVar.name());

      for (const auto & inputVarName : vvCovVars) {
        // Get input variable
        const oops::Variable inputVar = balVars_[inputVarName];

        // Get number of input levels
        const size_t nzJ = inputVar.getLevels();

        // Create vv-covariance field
        createField3D("vvCov", trans_->nw(), outputVar, inputVar, data_);

        // Get vv-covariance view
        auto vvCovView = getView3D("vvCov", outputVar, inputVar, data_);

        // Loop over ensemble members
        for (size_t je = 0; je < ne; ++je) {
          // Get member views
          const auto inputUnbalView = getView2D(inputVar, fsetEnsUnbal[je]);
          const auto outputUnbalView = getView2D(outputVar, fsetEnsUnbal[je]);

          // Update vv-covariance
          for (size_t js = 0; js < trans_->ns(); ++js) {
            for (size_t jw = 0; jw < trans_->nw(); ++jw) {
              if (trans_->includeWavenumber(js, jw)) {
                const double factor = trans_->spNorm(js);
                for (size_t jzI = 0; jzI < nzI; ++jzI) {
                  for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
                    vvCovView(jw, jzI, jzJ) += factor*inputUnbalView(js, jzJ)
                      *outputUnbalView(js, jzI);
                  }
                }
              }
            }
          }
        }

        // Get vv-covariance field
        auto vvCovField = getField("vvCov", outputVar, inputVar, data_);

        // Reduce and normalize vv-covariance
        trans_->reduceNormalizeCov(ne-1, vvCovField);

        // Filter vv-covariance
        trans_->filterCov(Lf_, vvCovField);
      }
    }
  }

  // Compute diagnostics
  computeDiagnostics();

  // Print norms
  print(oops::Log::test());

  oops::Log::trace() << classname() << "::directCalibration done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierBalance::iterativeCalibrationInit() {
  oops::Log::trace() << classname() << "::iterativeCalibrationInit starting" << std::endl;

  // Initialize iterative counters with zeroes
  iterativeN_ = 0;

  for (const auto & outputVar : balVarsExt_) {
    // Create perturbation field
    createField2D("pert", trans_->ns(), outputVar, data_);

    // Create mean field
    createField2D("mean", trans_->ns(), outputVar, data_);

    for (const auto & inputVar : xxCovVars(outputVar)) {
      // Create xx-covariance field
      createField3D("xxCov", trans_->nw(), outputVar, inputVar, data_);
    }
  }

  oops::Log::trace() << classname() << "::iterativeCalibrationInit done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierBalance::iterativeCalibrationUpdate(const oops::FieldSet3D & fset) {
  oops::Log::trace() << classname() << "::iterativeCalibrationUpdate starting" << std::endl;

  // Increment ensemble index (ie = ie + 1)
  iterativeN_++;

  // Sub-ensemble index
  const size_t ie = (params_.calibration.value()->subEnsSize.value() > 0) ?
    ((iterativeN_-1)%params_.calibration.value()->subEnsSize.value())+1 : iterativeN_;

  for (const auto & outputVar : balVarsExt_) {
    // Get number of output levels
    const size_t nzI = outputVar.getLevels();

    // Get member view
    const auto view = getView2D(outputVar, fset);

    // Get perturbation view
    auto outputPertView = getView2D("pert", outputVar, data_);

    // Get mean view
    auto meanView = getView2D("mean", outputVar, data_);

    if (ie == 1) {
      // Reset mean if a new sub-ensemble starts
      meanView.assign(0.0);
    }

    // Remove mean (pert = state - mean)
    for (size_t js = 0; js < trans_->ns(); ++js) {
      for (size_t jzI = 0; jzI < nzI; ++jzI) {
        outputPertView(js, jzI) = view(js, jzI) - meanView(js, jzI);
      }
    }

    for (const auto & inputVar : xxCovVars(outputVar)) {
      // Get number of input levels
      const size_t nzJ = inputVar.getLevels();

      // Get input perturbation view
      const auto inputPertView = getView2D("pert", inputVar, data_);

      if (ie > 1) {
        // Get xx-covariance view
        auto xxCovView = getView3D("xxCov", outputVar, inputVar, data_);

        // Update xx-covariance (cov = cov + (ie-1)/ie * inputPert * outputPert)
        for (size_t js = 0; js < trans_->ns(); ++js) {
          for (size_t jw = 0; jw < trans_->nw(); ++jw) {
            if (trans_->includeWavenumber(js, jw)) {
              const double factor = static_cast<double>(ie-1)/static_cast<double>(ie)
                *trans_->spNorm(js);
              for (size_t jzI = 0; jzI < nzI; ++jzI) {
                for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
                  xxCovView(jw, jzI, jzJ) += factor*inputPertView(js, jzJ)
                    *outputPertView(js, jzI);
                }
              }
            }
          }
        }
      }
    }

    // Update mean (mean = mean + 1 / ie * pert)
    const double factor = 1.0/static_cast<double>(ie);
    for (size_t js = 0; js < trans_->ns(); ++js) {
      for (size_t jzI = 0; jzI < nzI; ++jzI) {
        meanView(js, jzI) += factor*outputPertView(js, jzI);
      }
    }
  }

  oops::Log::trace() << classname() << "::iterativeCalibrationUpdate done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierBalance::iterativeCalibrationFinal() {
  oops::Log::trace() << classname() << "::iterativeCalibrationFinal starting" << std::endl;

  // Compute number of sub-ensembles
  size_t nSubEns = 1;
  if (params_.calibration.value()->subEnsSize.value() > 0) {
    // Check sub-ensembles size
    ASSERT(iterativeN_%params_.calibration.value()->subEnsSize.value() == 0);

    // Get number of sub-ensembles
    nSubEns = iterativeN_/params_.calibration.value()->subEnsSize.value();
  }

  for (const auto & outputVar : balVarsExt_) {
    for (const auto & inputVar : xxCovVars(outputVar)) {
      // Get xx-covariance field
      auto xxCovField = getField("xxCov", outputVar, inputVar, data_);

      // Reduce and normalize xx-covariance
      trans_->reduceNormalizeCov(iterativeN_-nSubEns, xxCovField);
    }
  }

  // Compute regressions from xx-covariances
  computeRegressionsFromCovariances();

  // Compute diagnostics
  computeDiagnostics();

  // Print norms
  print(oops::Log::test());

  oops::Log::trace() << classname() << "::iterativeCalibrationFinal done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierBalance::write() const {
  oops::Log::trace() << classname() << "::write starting" << std::endl;

  // Get number of variables
  size_t nVar = nCmp_;
  if (params_.write.value()) {
    if (params_.write.value()->writeCovariance.value()) {
      // Add covariance variables
      for (const auto & outputVar : balVarsExt_) {
        nVar += xxCovVars(outputVar).size();
      }
      for (const auto & row : params_.rows.value()) {
        if (row.inputVars.value().size() > 0) {
          ++nVar;
        }
      }
    }
    if (params_.calibration.value()) {
      // Add diagnostic variables
      for (const auto & row : params_.rows.value()) {
        nVar += 2;
        if (row.inputVars.value().size() > 0) {
          nVar += 2;
        }
        nVar += 4*row.inputVars.value().size();
      }
    }

    // NetCDF IDs
    int retval, ncId, nwGlbId;
    std::vector<int> nzIId(balVarsExt_.size()), nzJId(balVarsExt_.size());
    std::vector<int> dProfId(1), dVertCovId(2), dSpecId(2), dCovId(3), varId(nVar);

    // NetCDF file path
    const std::string ncFilePath = params_.write.value()->outputFile.value();

    // Definition mode
    size_t jVar = 0;

    if (comm_.rank() == 0) {
      // Create NetCDF file
      if ((retval = nc_create(ncFilePath.c_str(), NC_64BIT_OFFSET | NC_CLOBBER, &ncId)))
        ERR(retval, ncFilePath);

      // Create horizontal dimension
      if ((retval = nc_def_dim(ncId, "nwGlb", trans_->nwGlb(), &nwGlbId))) ERR(retval, "nwGlb");

      // Create vertical dimensions
      for (const auto & var : balVarsExt_) {
        // Get number of levels
        const size_t nz = var.getLevels();

        // Create dimensions
        const size_t jvar = balVarsExt_.find(var.name());
        const std::string nzIName = "nzI_" + var.name();
        const std::string nzJName = "nzJ_" + var.name();
        if ((retval = nc_def_dim(ncId, nzIName.c_str(), nz, &nzIId[jvar]))) ERR(retval, nzIName);
        if ((retval = nc_def_dim(ncId, nzJName.c_str(), nz, &nzJId[jvar]))) ERR(retval, nzJName);
      }

      // Dimensions arrays, horizontal part
      dSpecId[0] = nwGlbId;
      dCovId[0] = nwGlbId;

      for (const auto & row : params_.rows.value()) {
        // Get output variable
        const oops::Variable outputVar = balVarsExt_[row.outputVar.value()];

        for (const auto & inputVarName : row.inputVars.value()) {
          // Get input variable
          const oops::Variable inputVar = balVarsExt_[inputVarName];

          // Dimensions array, vertical part
          const size_t jvarI = balVarsExt_.find(outputVar.name());
          const size_t jvarJ = balVarsExt_.find(inputVar.name());
          dCovId[1] = nzIId[jvarI];
          dCovId[2] = nzJId[jvarJ];

          // Define variable
          defNcVar("reg", outputVar, inputVar, ncId, dCovId, varId, jVar);
        }
      }

      if (params_.write.value()->writeCovariance.value()) {
        // xx-covariances
        for (const auto & outputVar : balVarsExt_) {
          for (const auto & inputVar : xxCovVars(outputVar)) {
            // Dimensions array, vertical part
            const size_t jvarI = balVarsExt_.find(outputVar.name());
            const size_t jvarJ = balVarsExt_.find(inputVar.name());
            dCovId[1] = nzIId[jvarI];
            dCovId[2] = nzJId[jvarJ];

            // Define variable
            defNcVar("xxCov", outputVar, inputVar, ncId, dCovId, varId, jVar);
          }
        }

        // vv-covariances
        for (const auto & row : params_.rows.value()) {
          if (row.inputVars.value().size() > 0) {
            // Get variables
            const oops::Variable outputVar = balVarsExt_[row.outputVar.value()];

            // Dimensions array, vertical part
            const size_t jvarI = balVarsExt_.find(outputVar.name());
            const size_t jvarJ = jvarI;
            dCovId[1] = nzIId[jvarI];
            dCovId[2] = nzJId[jvarJ];

            // Define variable
            defNcVar("vvCov", outputVar, ncId, dCovId, varId, jVar);
          }
        }
      }

      if (params_.calibration.value()) {
        for (const auto & row : params_.rows.value()) {
          // Get output variable
          const oops::Variable outputVar = balVarsExt_[row.outputVar.value()];

          // Dimensions array, vertical part
          const size_t jvarI = balVarsExt_.find(outputVar.name());
          dProfId[0] = nzIId[jvarI];
          dVertCovId[0] = nzIId[jvarI];
          dSpecId[1] = nzIId[jvarI];

          // Define variables
          defNcVar("xxStdDevProf", outputVar, ncId, dProfId, varId, jVar);
          defNcVar("xxVarSpec", outputVar, ncId, dSpecId, varId, jVar);
          if (row.inputVars.value().size() > 0) {
            defNcVar("explVarRatioProf", outputVar, ncId, dProfId, varId, jVar);
            defNcVar("explVarRatioSpec", outputVar, ncId, dSpecId, varId, jVar);
          }

          for (const auto & inputVarName : row.inputVars.value()) {
            // Get input variable
            const oops::Variable inputVar = balVarsExt_[inputVarName];

            // Dimensions array, vertical part
            const size_t jvarJ = balVarsExt_.find(inputVar.name());
            dVertCovId[1] = nzJId[jvarJ];

            // Define variables
            const std::string varRatioProfName = "varRatioProf_from_" + inputVar.name() + "_for";
            defNcVar(varRatioProfName, outputVar, ncId, dProfId, varId, jVar);
            const std::string varRatioSpecName = "varRatioSpec_from_" + inputVar.name() + "_for";
            defNcVar(varRatioSpecName, outputVar, ncId, dSpecId, varId, jVar);
            defNcVar("xvVertCov", outputVar, inputVar, ncId, dVertCovId, varId, jVar);
            defNcVar("xvVertCor", outputVar, inputVar, ncId, dVertCovId, varId, jVar);
          }
        }
      }

      // End definition mode
      if ((retval = nc_enddef(ncId))) ERR(retval, ncFilePath);
    }

    // Data mode
    jVar = 0;

    for (const auto & row : params_.rows.value()) {
      // Get output variable
      const oops::Variable outputVar = balVarsExt_[row.outputVar.value()];

      for (const auto & inputVarName : row.inputVars.value()) {
        // Get input variable
        const oops::Variable inputVar = balVarsExt_[inputVarName];

        // Write variable
        putNcVar("reg", outputVar, inputVar, ncId, varId, jVar, data_, *trans_);
      }
    }

    if (params_.write.value()->writeCovariance.value()) {
      // xx-covariances
      for (const auto & outputVar : balVarsExt_) {
        for (const auto & inputVar : xxCovVars(outputVar)) {
          // Write variable
          putNcVar("xxCov", outputVar, inputVar, ncId, varId, jVar, data_, *trans_);
        }
      }

      // vv-covariances
      for (const auto & row : params_.rows.value()) {
        if (row.inputVars.value().size() > 0) {
          // Get variables
          const oops::Variable outputVar = balVarsExt_[row.outputVar.value()];

          // Write variable
          putNcVar("vvCov", outputVar, ncId, varId, jVar, data_, *trans_);
        }
      }
    }

    if (params_.calibration.value()) {
      for (const auto & row : params_.rows.value()) {
        // Get output variable
        const oops::Variable outputVar = balVarsExt_[row.outputVar.value()];

        // Write variables
        putNcVar("xxStdDevProf", outputVar, ncId, varId, jVar, data_, *trans_);
        putNcVar("xxVarSpec", outputVar, ncId, varId, jVar, data_, *trans_);
        if (row.inputVars.value().size() > 0) {
          putNcVar("explVarRatioProf", outputVar, ncId, varId, jVar, data_, *trans_);
          putNcVar("explVarRatioSpec", outputVar, ncId, varId, jVar, data_, *trans_);
        }

        for (const auto & inputVarName : row.inputVars.value()) {
          // Get input variable
          const oops::Variable inputVar = balVarsExt_[inputVarName];

          // Write variables
          const std::string varRatioProfName = "varRatioProf_from_" + inputVar.name() + "_for";
          putNcVar(varRatioProfName, outputVar, ncId, varId, jVar, data_, *trans_);
          const std::string varRatioSpecName = "varRatioSpec_from_" + inputVar.name() + "_for";
          putNcVar(varRatioSpecName, outputVar, ncId, varId, jVar, data_, *trans_);
          putNcVar("xvVertCov", outputVar, inputVar, ncId, varId, jVar, data_, *trans_);
          putNcVar("xvVertCor", outputVar, inputVar, ncId, varId, jVar, data_, *trans_);
        }
      }
    }

    if (comm_.rank() == 0) {
      // Close file
      if ((retval = nc_close(ncId))) ERR(retval, ncFilePath);
    }
  }

  oops::Log::trace() << classname() << "::write done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierBalance::readCovariance() {
  oops::Log::trace() << classname() << "::readCovariance starting" << std::endl;

  for (const auto & outputVar : balVarsExt_) {
    for (const auto & inputVar : xxCovVars(outputVar)) {
      // Create old xx-covariance field
      createField3D("xxOldCov", trans_->nw(), outputVar, inputVar, data_);
    }
  }

  // NetCDF file path
  const std::string ncFilePath = *params_.calibration.value()->oldCovInputFile.value();

  // NetCDF IDs
  int ncId, retval, nwGlbId, nzIId, nzJId, varId;
  size_t nwGlbFromFile, nzIFromFile, nzJFromFile;

  if (comm_.rank() == 0) {
    // Open NetCDF file
    if ((retval = nc_open(ncFilePath.c_str(), NC_NOWRITE, &ncId))) ERR(retval, ncFilePath);
  }

  // xx-covariances
  for (const auto & outputVar : balVarsExt_) {
    // Get number of output levels
    const size_t nzI = outputVar.getLevels();

    for (const auto & inputVar : xxCovVars(outputVar)) {
      // Get number of input levels
      const size_t nzJ = inputVar.getLevels();

      // Define global covariance vector
      std::vector<double> covVecGlb;

      if (comm_.rank() == 0) {
        // Check dimensions
        const std::string nzIName = "nzI_" + outputVar.name();
        const std::string nzJName = "nzJ_" + inputVar.name();
        if ((retval = nc_inq_dimid(ncId, "nwGlb", &nwGlbId))) ERR(retval, "nwGlb");
        if ((retval = nc_inq_dimid(ncId, nzIName.c_str(), &nzIId))) ERR(retval, nzIName);
        if ((retval = nc_inq_dimid(ncId, nzJName.c_str(), &nzJId))) ERR(retval, nzJName);
        if ((retval = nc_inq_dimlen(ncId, nwGlbId, &nwGlbFromFile))) ERR(retval, "nwGlb");
        if ((retval = nc_inq_dimlen(ncId, nzIId, &nzIFromFile))) ERR(retval, nzIName);
        if ((retval = nc_inq_dimlen(ncId, nzJId, &nzJFromFile))) ERR(retval, nzJName);
        ASSERT(nwGlbFromFile == trans_->nwGlb());
        ASSERT(nzIFromFile == nzI);
        ASSERT(nzJFromFile == nzJ);

        // Allocate global covariance vector
        covVecGlb.resize(trans_->nwGlb()*nzI*nzJ);

        // Read data
        const std::string xxCovFieldName = fieldName("xxCov", outputVar, inputVar);
        if ((retval = nc_inq_varid(ncId, xxCovFieldName.c_str(), &varId)))
          ERR(retval, xxCovFieldName);
        if ((retval = nc_get_var_double(ncId, varId, covVecGlb.data())))
          ERR(retval, xxCovFieldName);
      }

      // Get old xx-covariance field
      auto xxOldCovField = getField("xxOldCov", outputVar, inputVar, data_);

      // Scatter covariance vector
      trans_->scatterCov(covVecGlb, xxOldCovField);
    }
  }

  if (comm_.rank() == 0) {
    // Close file
    if ((retval = nc_close(ncId))) ERR(retval, ncFilePath);
  }

  oops::Log::trace() << classname() << "::readCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierBalance::computeRegression(const std::vector<std::string> & inputVars,
                                         const oops::Variable & outputVar) {
  oops::Log::trace() << classname() << "::computeRegression starting" << std::endl;

  if (inputVars.size() == 0) {
    // Noting to do
    oops::Log::trace() << classname() << "::computeRegression done" << std::endl;
    return;
  }

  // Initialize aggregated size
  size_t nzJagg = 0;

  for (const auto & inputVarName : inputVars) {
    // Get input variable
    const oops::Variable inputVar = balVars_[inputVarName];

    // Get number of levels
    const size_t nzJ = inputVar.getLevels();

    // Update aggregated number of levels
    nzJagg += nzJ;

    // Create regression field
    createField3D("reg", trans_->nw(), outputVar, inputVar, data_);
  }

  // Get number of levels
  const size_t nzI = outputVar.getLevels();

  // Create eigen matrices
  Eigen::MatrixXd vvCovMat(nzJagg, nzJagg);
  Eigen::MatrixXd xvCovMat(nzI, nzJagg);
  Eigen::MatrixXd regMat(nzI, nzJagg);

  for (size_t jw = 0; jw < trans_->nw(); ++jw) {
    // Get total wavenumber
    const size_t jwGlb = jw + trans_->nwStart();

    if (jwGlb == 0) {
      // No regression
      for (size_t jzI = 0; jzI < nzI; ++jzI) {
        for (size_t jzJ = 0; jzJ < nzJagg; ++jzJ) {
          regMat(jzI, jzJ) = 0.0;
        }
      }
    } else {
      // Initialize aggregated offset
      size_t dzJIagg = 0;

      for (const auto & inputVarIName : inputVars) {
        // Get input variable
        const oops::Variable inputVarI = balVars_[inputVarIName];

        // Get number of levels
        const size_t nzJI = inputVarI.getLevels();

        // Initialize aggregated offset
        size_t dzJJagg = 0;

        for (const auto & inputVarJName : inputVars) {
          if (balVars_.find(inputVarJName) <= balVars_.find(inputVarIName)) {
            // Get input variable
            const oops::Variable inputVarJ = balVars_[inputVarJName];

            // Get number of levels
            const size_t nzJJ = inputVarJ.getLevels();

            // Get vv-covariance view
            const auto vvCovView = getView3D("vvCov", inputVarI, inputVarJ, data_);

            // Convert to Eigen format
            for (size_t jzJI = 0; jzJI < nzJI; ++jzJI) {
              for (size_t jzJJ = 0; jzJJ < nzJJ; ++jzJJ) {
                vvCovMat(dzJIagg+jzJI, dzJJagg+jzJJ) = vvCovView(jw, jzJI, jzJJ);
                vvCovMat(dzJJagg+jzJJ, dzJIagg+jzJI) = vvCovMat(dzJIagg+jzJI, dzJJagg+jzJJ);
              }
            }

            // Update aggregated offset
            dzJJagg += nzJJ;
          }
        }

        // Get xv-covariance view
        const auto xvCovView = getView3D("xvCov", outputVar, inputVarI, data_);

        // Convert to Eigen format
        for (size_t jzI = 0; jzI < nzI; ++jzI) {
          for (size_t jzJ = 0; jzJ < nzJI; ++jzJ) {
            xvCovMat(jzI, dzJIagg+jzJ) = xvCovView(jw, jzI, jzJ);
          }
        }

        // Update aggregated offset
        dzJIagg += nzJI;
      }

      // Solve linear regression problem
      if (params_.calibration.value()->remainingVar.value() < 1.0) {
        // Split standard-deviation and correlation
        Eigen::VectorXd stdDevInv(nzJagg);
        for (size_t jzJ = 0; jzJ < nzJagg; ++jzJ) {
          ASSERT(vvCovMat(jzJ, jzJ) > 0.0);
          stdDevInv[jzJ] = 1.0/std::sqrt(vvCovMat(jzJ, jzJ));
        }
        for (size_t jzJI = 0; jzJI < nzJagg; ++jzJI) {
          for (size_t jzJJ = 0; jzJJ < nzJagg; ++jzJJ) {
            vvCovMat(jzJI, jzJJ) *= stdDevInv[jzJI]*stdDevInv[jzJJ];
          }
        }

        // Eigendecomposition to keep only a fraction of the spectrum variance
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
        es.compute(vvCovMat);
        Eigen::VectorXd d = es.eigenvalues();
        Eigen::MatrixXd V = es.eigenvectors();
        Eigen::MatrixXd vvCovInvMat = Eigen::MatrixXd::Zero(nzJagg, nzJagg);
        double spectrumVarianceTarget = 0.0;
        for (int jzJ = nzJagg-1; jzJ >= 0; --jzJ) {
          spectrumVarianceTarget += d[jzJ];
        }
        spectrumVarianceTarget *= params_.calibration.value()->remainingVar.value();
        double spectrumVariance = 0.0;
        for (int jzJ = nzJagg-1; jzJ >= 0; --jzJ) {
          spectrumVariance += d[jzJ];
          if (spectrumVariance > spectrumVarianceTarget) {
            d[jzJ] = 0.0;
          } else {
            ASSERT(d[jzJ] > 0.0);
            d[jzJ] = 1.0/d[jzJ];
          }
        }
        vvCovInvMat = stdDevInv.asDiagonal()*V*d.asDiagonal()*V.transpose()
          *stdDevInv.asDiagonal();
        regMat = xvCovMat*vvCovInvMat;
      } else {
        // Direct inversion
        regMat = (vvCovMat.selfadjointView<Eigen::Lower>().ldlt().solve(xvCovMat.transpose()))
          .transpose();
      }

      // Initialize aggregated offset
      size_t dzJagg = 0;

      for (const auto & inputVarName : inputVars) {
        // Get input variable
        const oops::Variable inputVar = balVars_[inputVarName];

        // Get number of levels
        const size_t nzJ = inputVar.getLevels();

        // Get regression view
        auto regView = getView3D("reg", outputVar, inputVar, data_);

        // Copy to fields
        for (size_t jzI = 0; jzI < nzI; ++jzI) {
          for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
            regView(jw, jzI, jzJ) = regMat(jzI, dzJagg+jzJ);
          }
        }

        // Update aggregated offset
        dzJagg += nzJ;
      }
    }
  }

  oops::Log::trace() << classname() << "::computeRegression done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierBalance::computeRegressionsFromCovariances() {
  oops::Log::trace() << classname() << "::computeRegressionsFromCovariances starting"
    << std::endl;

  // Combine xx-covariance with an old xx-covariance
  if (params_.calibration.value()) {
    if (params_.calibration.value()->oldCovInputFile.value()) {
      // Check parameters consistency
      ASSERT(params_.calibration.value()->halfLife.value());

      // Read old xx-covariance
      readCovariance();

      // Define update factor
      const double halfLife = *params_.calibration.value()->halfLife.value();
      const double alphaInf = 1.0-std::exp(-std::log(2.0)/halfLife);
      double updateFactor = alphaInf;
      if (params_.calibration.value()->cycleIndex.value()) {
        const size_t cycleIndex = *params_.calibration.value()->cycleIndex.value();
        ASSERT(cycleIndex > 0);
        updateFactor /= 1.0-std::pow(1.0-alphaInf, static_cast<double>(cycleIndex+1));
      }

      // xx-covariances
      for (const auto & outputVar : balVarsExt_) {
        // Get number of output levels
        const size_t nzI = outputVar.getLevels();

        for (const auto & inputVar : xxCovVars(outputVar)) {
          // Get number of input levels
          const size_t nzJ = inputVar.getLevels();

          // Get old xx-covariance view
          const auto xxOldCovView = getView3D("xxOldCov", outputVar, inputVar, data_);

          // Get xx-covariance view
          auto xxCovView = getView3D("xxCov", outputVar, inputVar, data_);

          // Combine xx-covariances
          for (size_t jw = 0; jw < trans_->nw(); ++jw) {
            for (size_t jzI = 0; jzI < nzI; ++jzI) {
              for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
                xxCovView(jw, jzI, jzJ) = updateFactor*xxCovView(jw, jzI, jzJ)
                  + (1.0-updateFactor)*xxOldCovView(jw, jzI, jzJ);
              }
            }
          }
        }
      }
    }
  }

  for (const auto & row : params_.rows.value()) {
    // Get output variable
    const oops::Variable outputVar = balVars_[row.outputVar.value()];

    // Get number of output levels
    const size_t nzI = outputVar.getLevels();

    for (const auto & inputVarName : row.inputVars.value()) {
      // Get input variable
      const oops::Variable inputVar = balVars_[inputVarName];

      // Get number of input levels
      const size_t nzJ = inputVar.getLevels();

      // Create xv-covariance field
      createField3D("xvCov", trans_->nw(), outputVar, inputVar, data_);

      // Get xv-covariance view
      auto xvCovView = getView3D("xvCov", outputVar, inputVar, data_);

      // Compute xv-covariance
      for (size_t jw = 0; jw < trans_->nw(); ++jw) {
        // Using documentation indices
        const size_t jvarI = balVars_.find(outputVar.name());
        const size_t jvarJ = balVars_.find(inputVar.name());

        for (size_t jvarK = 0; jvarK <= jvarJ; ++jvarK) {
          // Get number of levels
          const size_t nzK = balVars_[jvarK].getLevels();

          // Get xx-covariance view
          const auto xxCovView = getView3D("xxCov", balVars_[jvarI], balVars_[jvarK], data_);

          if (jvarK == jvarJ) {
            // Identity A matrix
            for (size_t jzI = 0; jzI < nzI; ++jzI) {
              for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
                xvCovView(jw, jzI, jzJ) += xxCovView(jw, jzI, jzJ);
              }
            }
          } else {
            // Get A matrix view
            const auto aView = getView3D("a", balVars_[jvarJ], balVars_[jvarK], data_);

            // Matrix multiplication
            for (size_t jzI = 0; jzI < nzI; ++jzI) {
              for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
                for (size_t jzK = 0; jzK < nzK; ++jzK) {
                  xvCovView(jw, jzI, jzJ) += xxCovView(jw, jzI, jzK)*aView(jw, jzJ, jzK);
                }
              }
            }
          }
        }
      }

      // Get xv-covariance field
      auto xvCovField = getField("xvCov", outputVar, inputVar, data_);

      // Filter xv-covariance
      trans_->filterCov(Lf_, xvCovField);
    }

    // Compute regression
    computeRegression(row.inputVars.value(), outputVar);

    // Compute A matrix
    for (const auto & inputVar : balVars_) {
      // Get number of input levels
      const size_t nzJ = inputVar.getLevels();

      if (balVars_.find(inputVar.name()) < balVars_.find(outputVar.name())) {
        // Create A matrix field
        createField3D("a", trans_->nw(), outputVar, inputVar, data_);

        // Get A matrix view
        auto aView = getView3D("a", outputVar, inputVar, data_);

        // Compute A matrix
        for (size_t jw = 0; jw < trans_->nw(); ++jw) {
          // Using documentation indices
          const size_t jvarI = balVars_.find(outputVar.name());
          const size_t jvarJ = balVars_.find(inputVar.name());
          for (size_t jvarK = jvarJ; jvarK <= jvarI-1; ++jvarK) {
            // Get number of levels
            const size_t nzK = balVars_[jvarK].getLevels();

            // Get regression field name
            const std::string regFieldName = fieldName("reg", outputVar, inputVar);

            if (data_.has(regFieldName)) {
              // Get regression view
              const auto regView = getView3D("reg", balVars_[jvarI], balVars_[jvarK], data_);

              if (jvarK == jvarJ) {
                // Identity temporary A matrix
                for (size_t jzI = 0; jzI < nzI; ++jzI) {
                  for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
                    aView(jw, jzI, jzJ) -= regView(jw, jzI, jzJ);
                  }
                }
              } else {
                // Get temporary A matrix view
                const auto tmpAView = getView3D("a", balVars_[jvarK], balVars_[jvarJ], data_);

                // Matrix multiplication
                for (size_t jzI = 0; jzI < nzI; ++jzI) {
                  for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
                    for (size_t jzK = 0; jzK < nzK; ++jzK) {
                      aView(jw, jzI, jzJ) -= regView(jw, jzI, jzK)*tmpAView(jw, jzK, jzJ);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    // vv-covariance variables
    std::vector<std::string> vvCovVars;
    if (outputVar.name() != balVars_.variables().back()) {
      vvCovVars = row.inputVars.value();
    }
    vvCovVars.push_back(outputVar.name());

    for (const auto & inputVarName : vvCovVars) {
      // Get input variable
      const oops::Variable inputVar = balVars_[inputVarName];

      // Get number of input levels
      const size_t nzJ = inputVar.getLevels();

      // Create vv-covariance field
      createField3D("vvCov", trans_->nw(), outputVar, inputVar, data_);

      // Get vv-covariance view
      auto vvCovView = getView3D("vvCov", outputVar, inputVar, data_);

      // Compute vv-covariance
      for (size_t jw = 0; jw < trans_->nw(); ++jw) {
        // Using documentation indices
        const size_t jvarI = balVars_.find(outputVar.name());
        const size_t jvarJ = balVars_.find(inputVarName);
        for (size_t jvarK = 0; jvarK <= jvarI; ++jvarK) {
          // Get number of levels
          const size_t nzK = balVars_[jvarK].getLevels();

          // Create temporary matrix field
          atlas::Field tmpField("tmp", make_datatype<double>(), make_shape(nzK, nzJ));

          // Get temporary matrix view
          auto tmpView = make_view<double, 2>(tmpField);

          // Initialize temporary matrix
          tmpView.assign(0.0);

          for (size_t jvarL = 0; jvarL <= jvarJ; ++jvarL) {
            // First product

            // Get number of levels
            const size_t nzL = balVars_[jvarL].getLevels();

            // Check whether the xx-covariance ajoint should be used
            const bool useAdjoint = jvarK < jvarL;

            // Get xx-covariance view
            const auto xxCovView = useAdjoint ?
              getView3D("xxCov", balVars_[jvarL], balVars_[jvarK], data_) :
              getView3D("xxCov", balVars_[jvarK], balVars_[jvarL], data_);

            // Update temporary matrix
            if (jvarL == jvarJ) {
              // Identity A matrix
              for (size_t jzK = 0; jzK < nzK; ++jzK) {
                for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
                  if (useAdjoint) {
                    tmpView(jzK, jzJ) += xxCovView(jw, jzJ, jzK);
                  } else {
                    tmpView(jzK, jzJ) += xxCovView(jw, jzK, jzJ);
                  }
                }
              }
            } else {
              // Get A matrix view
              const auto aView = getView3D("a", balVars_[jvarJ], balVars_[jvarL], data_);

              // Matrix multiplication
              for (size_t jzK = 0; jzK < nzK; ++jzK) {
                for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
                  for (size_t jzL = 0; jzL < nzL; ++jzL) {
                    if (useAdjoint) {
                      tmpView(jzK, jzJ) += xxCovView(jw, jzL, jzK)*aView(jw, jzJ, jzL);
                    } else {
                      tmpView(jzK, jzJ) += xxCovView(jw, jzK, jzL)*aView(jw, jzJ, jzL);
                    }
                  }
                }
              }
            }
          }

          // Second product
          if (jvarK == jvarI) {
            // Identity A matrix
            for (size_t jzI = 0; jzI < nzI; ++jzI) {
              for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
                vvCovView(jw, jzI, jzJ) += tmpView(jzI, jzJ);
              }
            }
          } else {
            // Get A matrix view
            const auto aView = getView3D("a", balVars_[jvarI], balVars_[jvarK], data_);

            // Matrix multiplication
            for (size_t jzI = 0; jzI < nzI; ++jzI) {
              for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
                for (size_t jzK = 0; jzK < nzK; ++jzK) {
                  vvCovView(jw, jzI, jzJ) += aView(jw, jzI, jzK)*tmpView(jzK, jzJ);
                }
              }
            }
          }
        }
      }

      // Get vv-covariance field
      auto vvCovField = getField("vvCov", outputVar, inputVar, data_);

      // Filter vv-covariance
      trans_->filterCov(Lf_, vvCovField);
    }
  }

  oops::Log::trace() << classname() << "::computeRegressionsFromCovariances done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierBalance::computeDiagnostics() {
  oops::Log::trace() << classname() << "::computeDiagnostics starting" << std::endl;

  // FieldSet for temporary covariances
  atlas::FieldSet fset;

  for (const auto & row : params_.rows.value()) {
    // Get output variable
    const oops::Variable outputVar = balVars_[row.outputVar.value()];

    // Get number of output levels
    const size_t nzI = outputVar.getLevels();

    // Create xx-covariance field
    createField3D("xxCov", trans_->nw(), outputVar, fset);

    // Get xx-covariance view
    auto xxCovView = getView3D("xxCov", outputVar, fset);

    for (const auto & inputVarName : row.inputVars.value()) {
      // Get input variable
      const oops::Variable inputVar = balVars_[inputVarName];

      // Get number of input levels
      const size_t nzJ = inputVar.getLevels();

      // Get regression view
      const auto regView = getView3D("reg", outputVar, inputVar, data_);

      // Get input vv-covariance view
      const auto inputVvCovView = getView3D("vvCov", inputVar, data_);

      // Balanced covariance name
      const std::string balCovName = "balCov_from_" + inputVar.name() + "_for";

      // Create balanced covariance field
      createField3D(balCovName, trans_->nw(), outputVar, fset);

      // Get balanced covariance view
      auto balCovView = getView3D(balCovName, outputVar, fset);

      // Compute balanced covariance
      for (size_t jw = 0; jw < trans_->nw(); ++jw) {
        for (size_t jzI1 = 0; jzI1 < nzI; ++jzI1) {
          for (size_t jzI2 = 0; jzI2 < nzI; ++jzI2) {
            for (size_t jzJ1 = 0; jzJ1 < nzJ; ++jzJ1) {
              for (size_t jzJ2 = 0; jzJ2 < nzJ; ++jzJ2) {
                balCovView(jw, jzI1, jzI2) += regView(jw, jzI1, jzJ1)*inputVvCovView(jw, jzJ1, jzJ2)
                  *regView(jw, jzI2, jzJ2);
              }
            }
          }
        }
      }

      // Add balanced covariance to xx-covariance
      for (size_t jw = 0; jw < trans_->nw(); ++jw) {
        for (size_t jzI1 = 0; jzI1 < nzI; ++jzI1) {
          for (size_t jzI2 = 0; jzI2 < nzI; ++jzI2) {
            xxCovView(jw, jzI1, jzI2) += balCovView(jw, jzI1, jzI2);
          }
        }
      }
    }

    // Get output vv-covariance view
    const auto outputVvCovView = getView3D("vvCov", outputVar, data_);

    // Add output vv-covariance to xx-covariance
    for (size_t jw = 0; jw < trans_->nw(); ++jw) {
      for (size_t jzI1 = 0; jzI1 < nzI; ++jzI1) {
        for (size_t jzI2 = 0; jzI2 < nzI; ++jzI2) {
          xxCovView(jw, jzI1, jzI2) += outputVvCovView(jw, jzI1, jzI2);
        }
      }
    }

    // Compute xx-variance profile
    std::vector<double> xxVarProf(nzI, 0.0);
    for (size_t js = 0; js < trans_->ns(); ++js) {
      if (trans_->q(js) == 0) {
        const double fact = trans_->spNorm(js)*trans_->spNorm(js);
        const size_t jw = trans_->jw(js);
        for (size_t jzI = 0; jzI < nzI; ++jzI) {
          xxVarProf[jzI] += xxCovView(jw, jzI, jzI)*fact;
        }
      }
    }

    // Reduce xx-variance profile
    comm_.allReduceInPlace(xxVarProf.begin(), xxVarProf.end(), eckit::mpi::sum());

    // Create xx-standard-deviation profile field
    createField1D("xxStdDevProf", outputVar, data_);

    // Get xx-standard-deviation profile view
    auto xxStdDevProfView = getView1D("xxStdDevProf", outputVar, data_);

    // Compute xx-standard deviation
    for (size_t jzI = 0; jzI < nzI; ++jzI) {
      xxStdDevProfView(jzI) = std::sqrt(xxVarProf[jzI]);
    }

    // Create xx-variance spectrum field
    createField2D("xxVarSpec", trans_->nw(), outputVar, data_);

    // xx-variance spectrum view
    auto xxVarSpecView = getView2D("xxVarSpec", outputVar, data_);

    // Compute xx-variance spectrum
    for (size_t jw = 0; jw < trans_->nw(); ++jw) {
      for (size_t jzI = 0; jzI < nzI; ++jzI) {
        xxVarSpecView(jw, jzI) = xxCovView(jw, jzI, jzI);
      }
    }

    for (const auto & inputVarName : row.inputVars.value()) {
      // Get input variable
      const oops::Variable inputVar = balVars_[inputVarName];

      // Get number of input levels
      const size_t nzJ = inputVar.getLevels();

      // Balanced covariance name
      const std::string balCovName = "balCov_from_" + inputVar.name() + "_for";

      // Get balanced covariance view
      const auto balCovView = getView3D(balCovName, outputVar, fset);

      // Compute balanced variance profile
      std::vector<double> balVarProf(nzI, 0.0);
      for (size_t js = 0; js < trans_->ns(); ++js) {
        if (trans_->q(js) == 0) {
          const double fact = trans_->spNorm(js)*trans_->spNorm(js);
          const size_t jw = trans_->jw(js);
          for (size_t jzI = 0; jzI < nzI; ++jzI) {
            balVarProf[jzI] += balCovView(jw, jzI, jzI)*fact;
          }
        }
      }

      // Reduce balanced variance profile
      comm_.allReduceInPlace(balVarProf.begin(), balVarProf.end(), eckit::mpi::sum());

      // Variance ratio profile name
      const std::string varRatioProfName = "varRatioProf_from_" + inputVar.name() + "_for";

      // Create variance ratio profile field
      createField1D(varRatioProfName, outputVar, data_);

      // Variance ratio profile view
      auto varRatioProfView = getView1D(varRatioProfName, outputVar, data_);

      // Compute variance ratio profile
      for (size_t jzI = 0; jzI < nzI; ++jzI) {
        if (xxVarProf[jzI] > 0.0) {
          varRatioProfView(jzI) = balVarProf[jzI]/xxVarProf[jzI];
        }
      }

      // Variance ratio spectrum name
      const std::string varRatioSpecName = "varRatioSpec_from_" + inputVar.name() + "_for";

      // Create variance ratio spectrum field
      createField2D(varRatioSpecName, trans_->nw(), outputVar, data_);

      // Variance ratio spectrum view
      auto varRatioSpecView = getView2D(varRatioSpecName, outputVar, data_);

      // Compute variance ratio spectrum
      for (size_t jw = 0; jw < trans_->nw(); ++jw) {
        for (size_t jzI = 0; jzI < nzI; ++jzI) {
          if (xxCovView(jw, jzI, jzI) > 0.0) {
            varRatioSpecView(jw, jzI) = balCovView(jw, jzI, jzI)/xxCovView(jw, jzI, jzI);
          }
        }
      }

      // Get xv-covariance view
      const auto xvCovView = getView3D("xvCov", outputVar, inputVar, data_);

      // Get input vv-covariance view
      const auto inputVvCovView = getView3D("vvCov", inputVar, data_);

      // Create vertical xv-covariance field
      createField2D("xvVertCov", outputVar, inputVar, data_);

      // Get vertical xv-covariance view
      auto xvVertCovView = getView2D("xvVertCov", outputVar, inputVar, data_);

      // Compute vertical xv-covariance
      for (size_t js = 0; js < trans_->ns(); ++js) {
        if (trans_->q(js) == 0) {
          const double fact = trans_->spNorm(js)*trans_->spNorm(js);
          const size_t jw = trans_->jw(js);
          for (size_t jzI = 0; jzI < nzI; ++jzI) {
            for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
              xvVertCovView(jzI, jzJ) += xvCovView(jw, jzI, jzJ)*fact;
            }
          }
        }
      }

      // Compute input vv-variance
      std::vector<double> vvVarProf(nzJ, 0.0);
      for (size_t js = 0; js < trans_->ns(); ++js) {
        if (trans_->q(js) == 0) {
          const double fact = trans_->spNorm(js)*trans_->spNorm(js);
          const size_t jw = trans_->jw(js);
          for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
            vvVarProf[jzJ] += inputVvCovView(jw, jzJ, jzJ)*fact;
          }
        }
      }

      // Reduce vv-variance profile
      comm_.allReduceInPlace(vvVarProf.begin(), vvVarProf.end(), eckit::mpi::sum());

      // Compute vv-standard-deviation
      std::vector<double> vvStdDevProf(nzJ, 0.0);
      for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
        vvStdDevProf[jzJ] = std::sqrt(vvVarProf[jzJ]);
      }

      // Create vertical xv-covariance and correlation fields
      createField2D("xvVertCor", outputVar, inputVar, data_);

      // Get vertical xv-covariance and correlation view
      auto xvVertCorView = getView2D("xvVertCor", outputVar, inputVar, data_);

      // Compute vertical xv-correlation
      for (size_t jzI = 0; jzI < nzI; ++jzI) {
        for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
          const double norm = xxStdDevProfView(jzI)*vvStdDevProf[jzJ];
          if (norm > 0.0) {
            xvVertCorView(jzI, jzJ) = xvVertCovView(jzI, jzJ)/norm;
          }
        }
      }
    }

    if (row.inputVars.value().size() > 0) {
      // Compute output vv-variance
      std::vector<double> vvVarProf(nzI, 0.0);
      for (size_t js = 0; js < trans_->ns(); ++js) {
        if (trans_->q(js) == 0) {
          const double fact = trans_->spNorm(js)*trans_->spNorm(js);
          const size_t jw = trans_->jw(js);
          for (size_t jzI = 0; jzI < nzI; ++jzI) {
            vvVarProf[jzI] += outputVvCovView(jw, jzI, jzI)*fact;
          }
        }
      }

      // Reduce vv-variance profile
      comm_.allReduceInPlace(vvVarProf.begin(), vvVarProf.end(), eckit::mpi::sum());

      // Compute vv-standard-deviation
      std::vector<double> vvStdDevProf(nzI, 0.0);
      for (size_t jzI = 0; jzI < nzI; ++jzI) {
        vvStdDevProf[jzI] = std::sqrt(vvStdDevProf[jzI]);
      }

      // Create explained variance ratio profile field
      createField1D("explVarRatioProf", outputVar, data_);

      // Get explained variance ratio profile view
      auto explVarRatioProfView = getView1D("explVarRatioProf", outputVar, data_);

      // Compute explained variance ratio profile
      for (size_t jzI = 0; jzI < nzI; ++jzI) {
        if (xxVarProf[jzI] > 0.0) {
          explVarRatioProfView(jzI) = 1.0 - vvVarProf[jzI]/xxVarProf[jzI];
        } else {
          explVarRatioProfView(jzI) = 0.0;
        }
      }

      // Create explained variance ratio spectrum field
      createField2D("explVarRatioSpec", trans_->nw(), outputVar, data_);

      // Explained variance ratio spectrum view
      auto explVarRatioSpecView = getView2D("explVarRatioSpec", outputVar, data_);

      // Compute explained variance ratio spectrum
      for (size_t jw = 0; jw < trans_->nw(); ++jw) {
        for (size_t jzI = 0; jzI < nzI; ++jzI) {
          if (xxCovView(jw, jzI, jzI) > 0.0) {
            explVarRatioSpecView(jw, jzI) = 1.0 - outputVvCovView(jw, jzI, jzI)
              /xxCovView(jw, jzI, jzI);
          }
        }
      }
    }
  }

  oops::Log::trace() << classname() << "::computeDiagnostics done" << std::endl;
}

// -----------------------------------------------------------------------------

oops::Variables BifourierBalance::xxCovVars(const oops::Variable & outputVar) const {
  oops::Log::trace() << classname() << "::xxCovVars starting" << std::endl;

  // Initialize input variables
  oops::Variables inputVars;

  if (balVars_.has(outputVar.name())) {
    // Normal variable: compute xx-covariances for all lower triangular blocks (including diagonal)
    for (const auto & inputVar : balVars_) {
      if (balVars_.find(inputVar.name()) <= balVars_.find(outputVar.name())) {
        inputVars.push_back(inputVar);
      }
    }
  } else {
    // Extra variable: compute xx-covariance with itself only
    inputVars.push_back(outputVar);
  }

  oops::Log::trace() << classname() << "::xxCovVars done" << std::endl;
  return inputVars;
}

// -----------------------------------------------------------------------------

void BifourierBalance::print(std::ostream & os) const {
  // Print norms
  os << "Regression norms: " << std::endl;
  for (const auto & row : params_.rows.value()) {
    // Get output variable
    const oops::Variable outputVar = balVars_[row.outputVar.value()];

    for (const auto & inputVarName : row.inputVars.value()) {
      // Get input variable
      const oops::Variable inputVar = balVars_[inputVarName];

      // Get regression field
      const auto regField = getField("reg", outputVar, inputVar, data_);

      // Print norms
      os << "- " << outputVar.name() << " from " << inputVar.name() << ": "
        << trans_->normCov(regField) << std::endl;
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
