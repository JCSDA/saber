/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/BifourierCovariance.h"

#include <Eigen/Dense>
#include <netcdf.h>

#include <algorithm>

#include "saber/bifourier/BifourierAromeLegacy.h"
#include "saber/bifourier/BifourierUtilities.h"

#define ERR(e, msg) {std::string s(nc_strerror(e)); \
  throw eckit::Exception(s + " : " + msg, Here());}

using atlas::array::make_datatype;
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
                                         const oops::FieldSet3D & fg) :
    SaberCentralBlockBase(params, xb.validTime()),
    comm_(gdata.comm()),
    activeVars_(activeVars),
    params_(params)
{
  oops::Log::trace() << classname() << "::BifourierCovariance starting" << std::endl;

  // Retrieve spectral transform
  trans_ = transStore_.retrieveTransform(gdata);

  // Create unified factor field
  atlas::Field unifiedFactorField("unifiedFactor", make_datatype<double>(),
    make_shape(trans_->ns()));

  // Get unified factor field view
  auto unifiedFactorView = make_view<double, 1>(unifiedFactorField);

  for (size_t js = 0; js < trans_->ns(); ++js) {
    // Unified factor for ellipses + FFT norm
    unifiedFactorView(js) = std::sqrt(trans_->spNorm(js)/trans_->normFFT());
  }

  // Add unified factor field
  data_.add(unifiedFactorField);

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
  trans_->createRandomSpectralFieldSet(fset.fieldSet(), activeVars_);

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

  // Get unified factor view
  const auto unifiedFactorView = make_view<double, 1>(data_["unifiedFactor"]);

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

    // Apply unified factor
    for (size_t js = 0; js < trans_->ns(); ++js) {
      for (size_t jz = 0; jz < nz; ++jz) {
        view(js, jz) *= unifiedFactorView(js);
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

  // Get unified factor view
  const auto unifiedFactorView = make_view<double, 1>(data_["unifiedFactor"]);

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

    // Apply unified factor
    for (size_t js = 0; js < trans_->ns(); ++js) {
      for (size_t jz = 0; jz < nz; ++jz) {
        copyView(js, jz) *= unifiedFactorView(js);
      }
    }

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

  // Read data
  if (params_.read.value()->inputFileFormat.value() == "arome legacy binary"
    || params_.read.value()->inputFileFormat.value() == "arome legacy netcdf") {
    arome_legacy::readCovariance(comm_, activeVars_, *params_.read.value(), *trans_, data_);
  } else {
    // NetCDF file path
    const std::string ncFilePath = params_.read.value()->inputFile.value();

    // NetCDF IDs
    int ncid, retval, nw_id, nzI_id, nzJ_id, corSqrt_id, stdDev_id;
    size_t nwGlbFromFile, nzIFromFile, nzJFromFile;

    if (comm_.rank() == 0) {
      // Open NetCDF file
      if ((retval = nc_open(ncFilePath.c_str(), NC_NOWRITE, &ncid))) ERR(retval, ncFilePath);
    }

    // Tag root
    size_t tagRoot = 123;

    for (const auto & var : activeVars_) {
      // Get number of levels
      const size_t nz = var.getLevels();

      // Get correlation square-root view
      auto corSqrtView = getView3D("corSqrt", var, data_);

      // Get standard-deviation view
      auto stdDevView = getViewProfile("stdDev", var, data_);

      // Define vector
      std::vector<double> corSqrtVec(trans_->nw()*nz*nz);
      std::vector<double> stdDevVec(nz);

      if (comm_.rank() == 0) {
        // Check dimensions
        const std::string nzIName = "nzI_" + var.name();
        const std::string nzJName = "nzJ_" + var.name();
        if ((retval = nc_inq_dimid(ncid, "nwGlb", &nw_id))) ERR(retval, "nwGlb");
        if ((retval = nc_inq_dimid(ncid, nzIName.c_str(), &nzI_id))) ERR(retval, nzIName);
        if ((retval = nc_inq_dimid(ncid, nzJName.c_str(), &nzJ_id))) ERR(retval, nzJName);
        if ((retval = nc_inq_dimlen(ncid, nw_id, &nwGlbFromFile))) ERR(retval, "nwGlb");
        if ((retval = nc_inq_dimlen(ncid, nzI_id, &nzIFromFile))) ERR(retval, nzIName);
        if ((retval = nc_inq_dimlen(ncid, nzJ_id, &nzJFromFile))) ERR(retval, nzJName);
        ASSERT(nwGlbFromFile == trans_->nwGlb());
        ASSERT(nzIFromFile == nz);
        ASSERT(nzJFromFile == nz);

        // Define global vector
        std::vector<double> corSqrtVecGlb(trans_->nwGlb()*nz*nz);

        // Get correlation square-root field name
        const std::string corSqrtFieldName = fieldName("corSqrt", var);

        // Get standard-deviation field name
        const std::string stdDevFieldName = fieldName("stdDev", var);

        // Get variables ID
        if ((retval = nc_inq_varid(ncid, corSqrtFieldName.c_str(), &corSqrt_id)))
          ERR(retval, corSqrtFieldName);
        if ((retval = nc_inq_varid(ncid, stdDevFieldName.c_str(), &stdDev_id)))
          ERR(retval, stdDevFieldName);

        // Read data
        if ((retval = nc_get_var_double(ncid, corSqrt_id, corSqrtVecGlb.data())))
          ERR(retval, corSqrtFieldName);
        if ((retval = nc_get_var_double(ncid, stdDev_id, stdDevVec.data())))
          ERR(retval, stdDevFieldName);

        // Copy data
        for (size_t jw = 0; jw < trans_->nw(); ++jw) {
          const size_t jwGlb = jw + trans_->nwStart();
          for (size_t jzI = 0; jzI < nz; ++jzI) {
            for (size_t jzJ = 0; jzJ < nz; ++jzJ) {
              const size_t index = jw*nz*nz + jzI*nz + jzJ;
              const size_t indexGlb = jwGlb*nz*nz + jzI*nz + jzJ;
              corSqrtVec[index] = corSqrtVecGlb[indexGlb];
            }
          }
        }

        // Send data
        int tag = tagRoot;
        for (size_t jt = 1; jt < comm_.size(); ++jt) {
          // Create vector
          std::vector<double> corSqrtVecSend(trans_->nwPerTask()[jt]*nz*nz);

          // Fill vector
          for (size_t jwSend = 0; jwSend < trans_->nwPerTask()[jt]; ++jwSend) {
            const size_t jwGlb = jwSend + trans_->nwStartPerTask()[jt];
            for (size_t jzI = 0; jzI < nz; ++jzI) {
              for (size_t jzJ = 0; jzJ < nz; ++jzJ) {
                const size_t indexSend = jwSend*nz*nz + jzI*nz + jzJ;
                const size_t indexGlb = jwGlb*nz*nz + jzI*nz + jzJ;
                corSqrtVecSend[indexSend] = corSqrtVecGlb[indexGlb];
              }
            }
          }

          // Send vector
          comm_.send(corSqrtVecSend.data(), trans_->nwPerTask()[jt]*nz*nz, jt, tag);

          // Update tag
          ++tag;
        }
      }

      if (comm_.rank() > 0) {
        // Receive vector
        comm_.receive(corSqrtVec.data(), trans_->nw()*nz*nz, 0, tagRoot+(comm_.rank()-1));
      }

      // MPI barrier
      comm_.barrier();

      // Broadcast data
      comm_.broadcast(stdDevVec.begin(), stdDevVec.end(), 0);

      // Deserialize
      for (size_t jw = 0; jw < trans_->nw(); ++jw) {
        for (size_t jzI = 0; jzI < nz; ++jzI) {
          for (size_t jzJ = 0; jzJ < nz; ++jzJ) {
            const size_t index = jw*nz*nz + jzI*nz + jzJ;
            corSqrtView(jw, jzI, jzJ) = corSqrtVec[index];
          }
        }
      }
      for (size_t jzJ = 0; jzJ < nz; ++jzJ) {
        stdDevView(jzJ) = stdDevVec[jzJ];
      }

      // Update tag root
      tagRoot += comm_.size();
    }

    if (comm_.rank() == 0) {
      // Close file
      if ((retval = nc_close(ncid))) ERR(retval, ncFilePath);
    }
  }

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierCovariance::write() const {
  oops::Log::trace() << classname() << "::write starting" << std::endl;

  // Check that output file is present
  ASSERT(params_.outputFile.value() != boost::none);

  // NetCDF IDs
  int retval, ncid, nw_id, nzI_id, nzJ_id, dStdDev_id[1], dCorSqrt_id[3],
    stdDev_id[activeVars_.size()], corSqrt_id[activeVars_.size()];

  // NetCDF file path
  const std::string ncFilePath = *params_.outputFile.value();

  // Definition mode
  size_t jvar = 0;

  if (comm_.rank() == 0) {
    // Create NetCDF file
    if ((retval = nc_create(ncFilePath.c_str(), NC_64BIT_OFFSET | NC_CLOBBER, &ncid)))
      ERR(retval, ncFilePath);

    // Create horizontal dimension
    if ((retval = nc_def_dim(ncid, "nwGlb", trans_->nwGlb(), &nw_id)))
      ERR(retval, "nwGlb");

    // Dimensions arrays, horizontal part
    dCorSqrt_id[0] = nw_id;

    for (const auto & var : activeVars_) {
      // Get number of levels
      const size_t nz = var.getLevels();

      // Create vertical dimensions
      const std::string nzIName = "nzI_" + var.name();
      const std::string nzJName = "nzJ_" + var.name();
      if ((retval = nc_def_dim(ncid, nzIName.c_str(), nz, &nzI_id))) ERR(retval, nzIName);
      if ((retval = nc_def_dim(ncid, nzJName.c_str(), nz, &nzJ_id))) ERR(retval, nzJName);

      // Dimensions array, vertical part
      dStdDev_id[0] = nzI_id;
      dCorSqrt_id[1] = nzI_id;
      dCorSqrt_id[2] = nzJ_id;

      // Get correlation square-root field name
      const std::string corSqrtFieldName = fieldName("corSqrt", var);

      // Get standard-deviation field name
      const std::string stdDevFieldName = fieldName("stdDev", var);

      // Define variables
      if ((retval = nc_def_var(ncid, corSqrtFieldName.c_str(), NC_DOUBLE, 3, dCorSqrt_id,
        &corSqrt_id[jvar]))) ERR(retval, corSqrtFieldName);
      if ((retval = nc_def_var(ncid, stdDevFieldName.c_str(), NC_DOUBLE, 1, dStdDev_id,
        &stdDev_id[jvar]))) ERR(retval, stdDevFieldName);

      // Update index
      ++jvar;
    }

    // End definition mode
    if ((retval = nc_enddef(ncid))) ERR(retval, ncFilePath);
  }

  // Data mode
  jvar = 0;

  for (const auto & var : activeVars_) {
    // Get number of levels
    const size_t nz = var.getLevels();

    // Get correlation square-root view
    const auto corSqrtView = getView3D("corSqrt", var, data_);

    // Get standard-deviation view
    const auto stdDevView = getViewProfile("stdDev", var, data_);

    // Define vectors
    std::vector<double> corSqrtVec(trans_->nwRoot()*nz*nz);
    std::vector<double> stdDevVec(nz);

    // Serialize
    for (size_t jw = 0; jw < trans_->nwRoot(); ++jw) {
      for (size_t jzI = 0; jzI < nz; ++jzI) {
        for (size_t jzJ = 0; jzJ < nz; ++jzJ) {
          const size_t index = jw*nz*nz + jzI*nz + jzJ;
          corSqrtVec[index] = corSqrtView(jw, jzI, jzJ);
        }
      }
    }
    for (size_t jz = 0; jz < nz; ++jz) {
      stdDevVec[jz] = stdDevView(jz);
    }

    // Define global vector
    std::vector<double> corSqrtVecGlb;
    if (comm_.rank() == 0) {
      corSqrtVecGlb.resize(trans_->nwGlb()*nz*nz);
    }

    // Gather vector
    std::vector<int> corSqrtCounts(trans_->wCounts());
    std::vector<int> corSqrtDispls(trans_->wDispls());
    for (size_t jt = 0; jt < comm_.size(); ++jt) {
      corSqrtCounts[jt] *= nz*nz;
      corSqrtDispls[jt] *= nz*nz;
    }
    comm_.gatherv(corSqrtVec.cbegin(), corSqrtVec.cend(), corSqrtVecGlb.begin(),
      corSqrtVecGlb.end(), corSqrtCounts, corSqrtDispls, 0);

    if (comm_.rank() == 0) {
      // Write data
      if ((retval = nc_put_var_double(ncid, corSqrt_id[jvar], corSqrtVecGlb.data())))
        ERR(retval, var.name() + "_corSqrt");
      if ((retval = nc_put_var_double(ncid, stdDev_id[jvar], stdDevVec.data())))
        ERR(retval, var.name() + "_stdDev");
      ++jvar;
    }
  }

  if (comm_.rank() == 0) {
    // Close file
    if ((retval = nc_close(ncid))) ERR(retval, ncFilePath);
  }

  oops::Log::trace() << classname() << "::write done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierCovariance::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber

