/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/BifourierAromeCovariance.h"

#include <netcdf.h>

#include <algorithm>
#include <vector>

#include "saber/bifourier/bifourier_arome_legacy.h"
#include "saber/bifourier/BifourierUtilities.h"

#define ERR(e, msg) {std::string s(nc_strerror(e)); \
  throw eckit::Exception(s + " : " + msg, Here());}

using atlas::array::make_datatype;
using atlas::array::make_shape;
using atlas::array::make_view;

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

static SaberCentralBlockMaker<BifourierAromeCovariance>
  makerBifourierAromeCovariance_("BifourierAromeCovariance");

// -----------------------------------------------------------------------------

void BifourierAromeCovariance::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;

  // Read data
  if (params_.read.value()->inputFileFormat.value() == "arome legacy binary"
    || params_.read.value()->inputFileFormat.value() == "arome legacy netcdf") {
    for (const auto & var : activeVars_) {
      // Create covariance field
      createField3D("cov", trans_->nw(), var, data_);
    }

    // Get number of levels
    const size_t nz = activeVars_["air_upward_absolute_vorticity"].getLevels();

    // Define global vectors
    std::vector<double> vorCovGlb;
    std::vector<double> divuCovGlb;
    std::vector<double> tPsuCovGlb;
    std::vector<double> quCovGlb;

    if (comm_.rank() == 0) {
      // Allocate global vectors
      vorCovGlb.resize(trans_->nwGlb()*nz*nz);
      divuCovGlb.resize(trans_->nwGlb()*nz*nz);
      tPsuCovGlb.resize(trans_->nwGlb()*(nz+1)*(nz+1));
      quCovGlb.resize(trans_->nwGlb()*nz*nz);

      if (params_.read.value()->inputFileFormat.value() == "arome legacy binary") {
        // Read Fortran unformatted file (from readjbdat96.F90)
        bifourier_arome_legacy_read_covariance_f90(params_.read.value()->toConfiguration(),
          trans_->nwGlb(), nz, vorCovGlb.data(), divuCovGlb.data(), tPsuCovGlb.data(),
           quCovGlb.data());
      } else if (params_.read.value()->inputFileFormat.value() == "arome legacy netcdf") {
        // NetCDF file path
        const std::string ncFilePath = params_.read.value()->inputFile.value();

        // NetCDF IDs
        int ncId, retval, dimId, varId;
        size_t nzFromFile, nwGlbFromFile;

        // Open NetCDF file
        if ((retval = nc_open(ncFilePath.c_str(), NC_NOWRITE, &ncId))) ERR(retval, ncFilePath);

        // Check dimensions
        if ((retval = nc_inq_dimid(ncId, "NFLEV", &dimId))) ERR(retval, "NFLEV");
        if ((retval = nc_inq_dimlen(ncId, dimId, &nzFromFile))) ERR(retval, "NFLEV");
        ASSERT(nzFromFile == nz);
        if ((retval = nc_inq_dimid(ncId, "NSMAXP1", &dimId))) ERR(retval, "NSMAXP1");
        if ((retval = nc_inq_dimlen(ncId, dimId, &nwGlbFromFile))) ERR(retval, "NSMAXP1");
        ASSERT(nwGlbFromFile == trans_->nwGlb());

        // Get variables
        if ((retval = nc_inq_varid(ncId, "VOR_VERTCOV", &varId))) ERR(retval, "VOR_VERTCOV");
        if ((retval = nc_get_var_double(ncId, varId, vorCovGlb.data()))) ERR(retval, "VOR_VERTCOV");
        if ((retval = nc_inq_varid(ncId, "DIVU_VERTCOV", &varId))) ERR(retval, "DIVU_VERTCOV");
        if ((retval = nc_get_var_double(ncId, varId, divuCovGlb.data())))
          ERR(retval, "DIVU_VERTCOV");
        if ((retval = nc_inq_varid(ncId, "TPSU_VERTCOV", &varId))) ERR(retval, "TPSU_VERTCOV");
        if ((retval = nc_get_var_double(ncId, varId, tPsuCovGlb.data())))
          ERR(retval, "TPSU_VERTCOV");
        if ((retval = nc_inq_varid(ncId, "QU_VERTCOV", &varId))) ERR(retval, "QU_VERTCOV");
        if ((retval = nc_get_var_double(ncId, varId, quCovGlb.data()))) ERR(retval, "QU_VERTCOV");

        // Close file
        if ((retval = nc_close(ncId))) ERR(retval, ncFilePath);
      }
    }

    // Scatter data
    for (const auto & var : activeVars_) {
      // Get covariance field
      auto covField = getField("cov", var, data_);

      // Scatter global vector
      if (var.name() == "air_upward_absolute_vorticity") {
        trans_->scatterCov(vorCovGlb, covField, true);
      }
      if (var.name() == "air_horizontal_divergence") {
        trans_->scatterCov(divuCovGlb, covField, true);
      }
      if (var.name() == "air_temperature_and_log_of_air_pressure_at_surface") {
        trans_->scatterCov(tPsuCovGlb, covField, true);
      }
      if (var.name() == "water_vapor_mixing_ratio_wrt_moist_air") {
        trans_->scatterCov(quCovGlb, covField, true);
      }
    }

    // Rescale covariance from AROME to block standard
    for (const auto & var : activeVars_) {
      // Get number of levels
      const size_t nz = var.getLevels();

      // Get covariance view
      auto covView = getView3D("cov", var, data_);

      for (size_t jw = 0; jw < trans_->nw(); ++jw) {
        // Get AROME weight
        const double zWeight = 1.0/aromeWeight(jw);

        // Apply weight
        for (size_t jzI = 0; jzI < nz; ++jzI) {
          for (size_t jzJ = 0; jzJ < nz; ++jzJ) {
            covView(jw, jzI, jzJ) *= zWeight;
          }
        }
      }
    }

    // Compute square-root
    computeSquareRoot();

    // Print norms
    print(oops::Log::test());
  } else {
    // Generic reader
    BifourierCovariance::read();
  }

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierAromeCovariance::write() const {
  oops::Log::trace() << classname() << "::write starting" << std::endl;

  if (params_.write.value() != boost::none) {
    // Write data
    if (params_.write.value()->outputFileFormat.value() == "arome legacy binary"
      || params_.write.value()->outputFileFormat.value() == "arome legacy netcdf") {
      // Create AROME covariance fieldset
      atlas::FieldSet aromeCovData;

      // Compute covariance from correlation square-root and standard-deviation if it is missing
      computeCovariance(aromeCovData);

      // Define global vectors
      std::vector<double> vorCovGlb;
      std::vector<double> divuCovGlb;
      std::vector<double> tPsuCovGlb;
      std::vector<double> quCovGlb;


      for (const auto & var : activeVars_) {
        // Get number of levels
        const size_t nz = var.getLevels();

        // Get covariance view
        const auto covView = getView3D("cov", var, aromeCovData);

        // Create AROME covariance field
        createField3D("aromeCov", trans_->nw(), var, aromeCovData);

        // Get AROME covariance view
        auto aromeCovView = getView3D("aromeCov", var, aromeCovData);

        for (size_t jw = 0; jw < trans_->nw(); ++jw) {
          // Get AROME weight
          const double zWeight = aromeWeight(jw);

          // Apply weight
          for (size_t jzI = 0; jzI < nz; ++jzI) {
            for (size_t jzJ = 0; jzJ < nz; ++jzJ) {
              aromeCovView(jw, jzI, jzJ) = covView(jw, jzI, jzJ)*zWeight;
            }
          }
        }
      }

      // Write AROME covariances
      for (const auto & var : activeVars_) {
        // Get covariance field
        const auto aromeCovField = getField("aromeCov", var, aromeCovData);

        // Gather covariance vector
        if (var.name() == "air_upward_absolute_vorticity") {
          trans_->gatherCov(aromeCovField, vorCovGlb, true);
        }
        if (var.name() == "air_horizontal_divergence") {
          trans_->gatherCov(aromeCovField, divuCovGlb, true);
        }
        if (var.name() == "air_temperature_and_log_of_air_pressure_at_surface") {
          trans_->gatherCov(aromeCovField, tPsuCovGlb, true);
        }
        if (var.name() == "water_vapor_mixing_ratio_wrt_moist_air") {
          trans_->gatherCov(aromeCovField, quCovGlb, true);
        }
      }

      if (comm_.rank() == 0) {
        // Get number of levels
        const size_t nz = activeVars_["air_upward_absolute_vorticity"].getLevels();

        if (params_.write.value()->outputFileFormat.value() == "arome legacy binary") {
          // Write Fortran unformatted file (from ewgsacov.F90)
          bifourier_arome_legacy_write_covariance_f90(params_.write.value()->toConfiguration(),
            trans_->nwGlb(), nz, vorCovGlb.data(), divuCovGlb.data(), tPsuCovGlb.data(),
             quCovGlb.data());
        } else if (params_.write.value()->outputFileFormat.value() == "arome legacy netcdf") {
          // NetCDF file path
          const std::string ncFilePath = params_.write.value()->outputFile.value();

          // NetCDF IDs
          int ncId, retval, nzId, nzP1Id, nwGlbId, dNzId[3], dNzP1Id[3],
            vorCovId, divuCovId, tPsuCovId, quCovId;

          // Create NetCDF file
          if ((retval = nc_create(ncFilePath.c_str(), NC_64BIT_OFFSET | NC_CLOBBER, &ncId)))
            ERR(retval, ncFilePath);

          // Create dimensions
          if ((retval = nc_def_dim(ncId, "NFLEV", nz, &nzId))) ERR(retval, "NFLEV");
          if ((retval = nc_def_dim(ncId, "NFLEVP1", nz+1, &nzP1Id))) ERR(retval, "NFLEVP1");
          if ((retval = nc_def_dim(ncId, "NSMAXP1", trans_->nwGlb(), &nwGlbId)))
            ERR(retval, "NSMAXP1");

          // Dimensions arrays
          dNzId[0] = nwGlbId;
          dNzId[1] = nzId;
          dNzId[2] = nzId;
          dNzP1Id[0] = nwGlbId;
          dNzP1Id[1] = nzP1Id;
          dNzP1Id[2] = nzP1Id;

          // Create variables
          if ((retval = nc_def_var(ncId, "VOR_VERTCOV", NC_DOUBLE, 3, dNzId, &vorCovId)))
            ERR(retval, "VOR_VERTCOV");
          if ((retval = nc_def_var(ncId, "DIVU_VERTCOV", NC_DOUBLE, 3, dNzId, &divuCovId)))
            ERR(retval, "DIVU_VERTCOV");
          if ((retval = nc_def_var(ncId, "TPSU_VERTCOV", NC_DOUBLE, 3, dNzP1Id, &tPsuCovId)))
            ERR(retval, "TPSU_VERTCOV");
          if ((retval = nc_def_var(ncId, "QU_VERTCOV", NC_DOUBLE, 3, dNzId, &quCovId)))
            ERR(retval, "QU_VERTCOV");

          // End definition mode
          if ((retval = nc_enddef(ncId))) ERR(retval, ncFilePath);

          // Write data
          if ((retval = nc_put_var_double(ncId, vorCovId, vorCovGlb.data())))
            ERR(retval, "VOR_VERTCOV");
          if ((retval = nc_put_var_double(ncId, divuCovId, divuCovGlb.data())))
            ERR(retval, "DIVU_VERTCOV");
          if ((retval = nc_put_var_double(ncId, tPsuCovId, tPsuCovGlb.data())))
            ERR(retval, "TPSU_VERTCOV");
          if ((retval = nc_put_var_double(ncId, quCovId, quCovGlb.data())))
            ERR(retval, "QU_VERTCOV");

          // Close file
          if ((retval = nc_close(ncId))) ERR(retval, ncFilePath);
        }
      }
    } else {
      // Generic writer
      BifourierCovariance::write();
    }
  }

  oops::Log::trace() << classname() << "::write done" << std::endl;
}

// -----------------------------------------------------------------------------

double BifourierAromeCovariance::aromeWeight(const size_t & jw) const {
  oops::Log::trace() << classname() << "::aromeWeight starting" << std::endl;

  // Get global total wavenumber
  const size_t jwGlb = jw + trans_->nwStart();

  // Constant coefficient
  const double zmovern = static_cast<double>(trans_->ellips().size())
    / static_cast<double>(trans_->nwGlb()-1);

  // Compute weight
  double zWeight;
  if (jwGlb != 0 && jwGlb != trans_->nwGlb()-1) {
    zWeight = 2.0*M_PI*static_cast<double>(jwGlb)*zmovern;
  } else if (jwGlb == 0) {
//    zWeight = M_PI*zmovern/4.0;
    zWeight = M_PI*zmovern/2.0;
  } else if (jwGlb == trans_->nwGlb()-1) {
    zWeight = M_PI*(static_cast<double>(trans_->nwGlb()-1)-0.25)*zmovern;
  }

  // REDNMC factor
  zWeight /= params_.rednmc.value()*params_.rednmc.value();

  oops::Log::trace() << classname() << "::aromeWeight starting" << std::endl;
  return zWeight;
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber

