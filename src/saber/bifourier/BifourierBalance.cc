/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/BifourierBalance.h"

#include <Eigen/Dense>
#include <netcdf.h>

#include <algorithm>

#include "oops/util/FieldSetOperations.h"

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

static SaberOuterBlockMaker<BifourierBalance> makerBifourierBalance_("BifourierBalance");

// -----------------------------------------------------------------------------

BifourierBalance::BifourierBalance(const oops::GeometryData & outerGeometryData,
                                   const oops::Variables & outerVars,
                                   const eckit::Configuration & covarConfig,
                                   const Parameters_ & params,
                                   const oops::FieldSet3D & xb,
                                   const oops::FieldSet3D & fg)
  : SaberOuterBlockBase(params, xb.validTime()),
    innerGeometryData_(outerGeometryData),
    comm_(outerGeometryData.comm()),
    innerVars_(outerVars),
    params_(params)
{
  oops::Log::trace() << classname() << "::BifourierBalance starting" << std::endl;

  // Retrieve spectral transform
  trans_ = transStore_.retrieveTransform(outerGeometryData);

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

  // Read data
  if (params_.read.value()->inputFileFormat.value() == "arome legacy binary"
    || params_.read.value()->inputFileFormat.value() == "arome legacy netcdf") {
    arome_legacy::readBalance(comm_, balVars_, *params_.read.value(), *trans_, data_);
  } else {
    // NetCDF file path
    const std::string ncFilePath = params_.read.value()->inputFile.value();

    // NetCDF IDs
    int ncid, retval, nwGlb_id, nzI_id, nzJ_id, varid;
    size_t nwGlbFromFile, nzIFromFile, nzJFromFile;

    if (comm_.rank() == 0) {
      // Open NetCDF file
      if ((retval = nc_open(ncFilePath.c_str(), NC_NOWRITE, &ncid))) ERR(retval, ncFilePath);
    }

    // Tag root
    size_t tagRoot = 123;

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

        // Get regression view
        auto regView = getView3D("reg", outputVar, inputVar, data_);

        // Define vector
        std::vector<double> vec(trans_->nw()*nzI*nzJ);

        if (comm_.rank() == 0) {
          // Check dimensions
          const std::string nzIName = "nzI_" + outputVar.name();
          const std::string nzJName = "nzJ_" + inputVar.name();
          if ((retval = nc_inq_dimid(ncid, "nwGlb", &nwGlb_id))) ERR(retval, "nwGlb");
          if ((retval = nc_inq_dimid(ncid, nzIName.c_str(), &nzI_id))) ERR(retval, nzIName);
          if ((retval = nc_inq_dimid(ncid, nzJName.c_str(), &nzJ_id))) ERR(retval, nzJName);
          if ((retval = nc_inq_dimlen(ncid, nwGlb_id, &nwGlbFromFile))) ERR(retval, "nwGlb");
          if ((retval = nc_inq_dimlen(ncid, nzI_id, &nzIFromFile))) ERR(retval, nzIName);
          if ((retval = nc_inq_dimlen(ncid, nzJ_id, &nzJFromFile))) ERR(retval, nzJName);
          ASSERT(nwGlbFromFile == trans_->nwGlb());
          ASSERT(nzIFromFile == nzI);
          ASSERT(nzJFromFile == nzJ);

          // Define global vector
          std::vector<double> vecGlb(trans_->nwGlb()*nzI*nzJ);

          // Read data
          const std::string regFieldName = fieldName("reg", outputVar, inputVar);
          if ((retval = nc_inq_varid(ncid, regFieldName.c_str(), &varid)))
            ERR(retval, regFieldName);
          if ((retval = nc_get_var_double(ncid, varid, vecGlb.data()))) ERR(retval, regFieldName);

          // Copy data
          for (size_t jw = 0; jw < trans_->nw(); ++jw) {
            const size_t jwGlb = jw + trans_->nwStart();
            for (size_t jzI = 0; jzI < nzI; ++jzI) {
              for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
                const size_t jj = jw*nzI*nzJ + jzI*nzJ + jzJ;
                const size_t jjGlb = jwGlb*nzI*nzJ + jzI*nzJ + jzJ;
                vec[jj] = vecGlb[jjGlb];
              }
            }
          }

          // Send data
          int tag = tagRoot;
          for (size_t jt = 1; jt < comm_.size(); ++jt) {
            // Create vector
            std::vector<double> vecSend(trans_->nwPerTask()[jt]*nzI*nzJ);

            // Fill vector
            for (size_t jwSend = 0; jwSend < trans_->nwPerTask()[jt]; ++jwSend) {
              const size_t jwGlb = jwSend + trans_->nwStartPerTask()[jt];
              for (size_t jzI = 0; jzI < nzI; ++jzI) {
                for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
                  const size_t jjSend = jwSend*nzI*nzJ + jzI*nzJ + jzJ;
                  const size_t jjGlb = jwGlb*nzI*nzJ + jzI*nzJ + jzJ;
                  vecSend[jjSend] = vecGlb[jjGlb];
                }
              }
            }

            // Send vector
            comm_.send(vecSend.data(), trans_->nwPerTask()[jt]*nzI*nzJ, jt, tag);

            // Update tag
            ++tag;
          }
        }

        if (comm_.rank() > 0) {
          // Receive vector
          comm_.receive(vec.data(), trans_->nw()*nzI*nzJ, 0, tagRoot+(comm_.rank()-1));
        }

        // MPI barrier
        comm_.barrier();

        // Deserialize
        for (size_t jw = 0; jw < trans_->nw(); ++jw) {
          for (size_t jzI = 0; jzI < nzI; ++jzI) {
            for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
              const size_t jj = jw*nzI*nzJ + jzI*nzJ + jzJ;
              regView(jw, jzI, jzJ) = vec[jj];
            }
          }
        }

        // Update tag root
        tagRoot += comm_.size();
      }
    }

    if (comm_.rank() == 0) {
      // Close file
      if ((retval = nc_close(ncid))) ERR(retval, ncFilePath);
    }
  }

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierBalance::write() const {
  oops::Log::trace() << classname() << "::write starting" << std::endl;

  // Check that output file is present
  ASSERT(params_.outputFile.value() != boost::none);

  // NetCDF IDs
  int retval, ncid, nwGlb_id, d3D_id[3], nzI_id[balVars_.size()], nzJ_id[balVars_.size()],
    reg_id[nCmp_];

  // NetCDF file path
  const std::string ncFilePath = *params_.outputFile.value();

  // Definition mode
  size_t jCmp = 0;

  if (comm_.rank() == 0) {
    // Create NetCDF file
    if ((retval = nc_create(ncFilePath.c_str(), NC_64BIT_OFFSET | NC_CLOBBER, &ncid)))
      ERR(retval, ncFilePath);

    // Create horizontal dimension
    if ((retval = nc_def_dim(ncid, "nwGlb", trans_->nwGlb(), &nwGlb_id))) ERR(retval, "nwGlb");

    // Create vertical dimensions
    for (const auto & var : balVars_) {
      // Get number of levels
      const size_t nz = var.getLevels();

      // Create dimensions
      const size_t jvar = balVars_.find(var.name());
      const std::string nzIName = "nzI_" + var.name();
      const std::string nzJName = "nzJ_" + var.name();
      if ((retval = nc_def_dim(ncid, nzIName.c_str(), nz, &nzI_id[jvar]))) ERR(retval, nzIName);
      if ((retval = nc_def_dim(ncid, nzJName.c_str(), nz, &nzJ_id[jvar]))) ERR(retval, nzJName);
    }

    // Dimensions arrays, horizontal part
    d3D_id[0] = nwGlb_id;

    for (const auto & row : params_.rows.value()) {
      // Get output variable
      const oops::Variable outputVar = balVars_[row.outputVar.value()];

      for (const auto & inputVarName : row.inputVars.value()) {
        // Get input variable
        const oops::Variable inputVar = balVars_[inputVarName];

        // Dimensions array, vertical part
        const size_t jvarI = balVars_.find(outputVar.name());
        const size_t jvarJ = balVars_.find(inputVar.name());
        d3D_id[1] = nzI_id[jvarI];
        d3D_id[2] = nzJ_id[jvarJ];

        // Define variable
        const std::string regFieldName = fieldName("reg", outputVar, inputVar);
        if ((retval = nc_def_var(ncid, regFieldName.c_str(), NC_DOUBLE, 3, d3D_id, &reg_id[jCmp])))
          ERR(retval, regFieldName);
        ++jCmp;
      }
    }

    // End definition mode
    if ((retval = nc_enddef(ncid))) ERR(retval, ncFilePath);
  }

  // Data mode
  jCmp = 0;

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

      // Get regression view
      auto regView = getView3D("reg", outputVar, inputVar, data_);

      // Serialize
      std::vector<double> vec(trans_->nwRoot()*nzI*nzJ);
      for (size_t jw = 0; jw < trans_->nwRoot(); ++jw) {
        for (size_t jzI = 0; jzI < nzI; ++jzI) {
          for (size_t jzJ = 0; jzJ < nzJ; ++jzJ) {
            const size_t jj = jw*nzI*nzJ + jzI*nzJ + jzJ;
            vec[jj] = regView(jw, jzI, jzJ);
          }
        }
      }

      // Define global vector
      std::vector<double> vecGlb;
      if (comm_.rank() == 0) {
        vecGlb.resize(trans_->nwGlb()*nzI*nzJ);
      }

      // Gather vector
      std::vector<int> wCounts(trans_->wCounts());
      std::vector<int> wDispls(trans_->wDispls());
      for (size_t jt = 0; jt < comm_.size(); ++jt) {
        wCounts[jt] *= nzI*nzJ;
        wDispls[jt] *= nzI*nzJ;
      }
      comm_.gatherv(vec.cbegin(), vec.cend(), vecGlb.begin(), vecGlb.end(), wCounts, wDispls, 0);

      if (comm_.rank() == 0) {
        // Write data
        const std::string regFieldName = fieldName("reg", outputVar, inputVar);
        if ((retval = nc_put_var_double(ncid, reg_id[jCmp], vecGlb.data())))
          ERR(retval, regFieldName);
        ++jCmp;
      }
    }
  }

  if (comm_.rank() == 0) {
    // Close file
    if ((retval = nc_close(ncid))) ERR(retval, ncFilePath);
  }

  oops::Log::trace() << classname() << "::write done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierBalance::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
