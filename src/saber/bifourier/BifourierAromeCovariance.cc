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

static CentralBlockMaker<BifourierAromeCovariance>
  makerBifourierAromeCovariance_("BifourierAromeCovariance");

// -----------------------------------------------------------------------------

void BifourierAromeCovariance::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;

  // Read data
  if (params_.read.value()->inputFileFormat.value() == "arome legacy binary"
    || params_.read.value()->inputFileFormat.value() == "arome legacy netcdf") {
    for (const auto & var : centralVars()) {
      // Create covariance field
      createField3D("cov", covar_->trans()->nw(), var, covar_->data());
    }

    // Get number of levels
    const size_t nz = centralVars()["air_upward_absolute_vorticity"].getLevels();

    // Define global vectors
    std::vector<double> vorCovGlb;
    std::vector<double> divuCovGlb;
    std::vector<double> tPsuCovGlb;
    std::vector<double> quCovGlb;

    // Define attributes
    eckit::LocalConfiguration attributes;
    attributes.set("nsmax", covar_->trans()->nwGlb()-1);
    attributes.set("nflev", nz);

    if (covar_->comm().rank() == 0) {
      // Allocate global vectors
      vorCovGlb.resize(covar_->trans()->nwGlb()*nz*nz);
      divuCovGlb.resize(covar_->trans()->nwGlb()*nz*nz);
      tPsuCovGlb.resize(covar_->trans()->nwGlb()*(nz+1)*(nz+1));
      quCovGlb.resize(covar_->trans()->nwGlb()*nz*nz);

      if (params_.read.value()->inputFileFormat.value() == "arome legacy binary") {
        // Read Fortran unformatted file (from readjbdat96.F90)
        const int nsmax = attributes.getDouble("nsmax");
        const int nflev = attributes.getDouble("nflev");
        bifourier_arome_legacy_read_covariance_f90(params_.read.value()->toConfiguration(),
          attributes, nsmax+1, nflev, vorCovGlb.data(), divuCovGlb.data(),
          tPsuCovGlb.data(), quCovGlb.data());
      } else if (params_.read.value()->inputFileFormat.value() == "arome legacy netcdf") {
        // NetCDF file path
        const std::string ncFilePath = params_.read.value()->inputFile.value();

        // NetCDF IDs
        int ncId, retval, dimId, varId;
        size_t nflevFile, nsmaxp1File;

        // Open NetCDF file
        if ((retval = nc_open(ncFilePath.c_str(), NC_NOWRITE, &ncId))) ERR(retval, ncFilePath);

        // Check dimensions
        if ((retval = nc_inq_dimid(ncId, "NFLEV", &dimId))) ERR(retval, "NFLEV");
        if ((retval = nc_inq_dimlen(ncId, dimId, &nflevFile))) ERR(retval, "NFLEV");
        ASSERT(nflevFile == nz);
        if ((retval = nc_inq_dimid(ncId, "NSMAXP1", &dimId))) ERR(retval, "NSMAXP1");
        if ((retval = nc_inq_dimlen(ncId, dimId, &nsmaxp1File))) ERR(retval, "NSMAXP1");
        ASSERT(nsmaxp1File == covar_->trans()->nwGlb());

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
    for (const auto & var : centralVars()) {
      // Get covariance field
      auto covField = getField("cov", var, covar_->data());

      // Scatter global vector
      if (var.name() == "air_upward_absolute_vorticity") {
        covar_->trans()->scatterCov(vorCovGlb, covField, true);
      }
      if (var.name() == "air_horizontal_divergence") {
        covar_->trans()->scatterCov(divuCovGlb, covField, true);
      }
      if (var.name() == "air_temperature_and_log_of_air_pressure_at_surface") {
        covar_->trans()->scatterCov(tPsuCovGlb, covField, true);
      }
      if (var.name() == "water_vapor_mixing_ratio_wrt_moist_air") {
        covar_->trans()->scatterCov(quCovGlb, covField, true);
      }
    }

    // Rescale covariance from AROME to block standard
    for (const auto & var : centralVars()) {
      // Get number of levels
      const size_t nz = var.getLevels();

      // Get covariance view
      auto covView = getView3D("cov", var, covar_->data());

      for (size_t jw = 0; jw < covar_->trans()->nw(); ++jw) {
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
    covar_->computeSquareRoot();

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

  if (params_.write.value()) {
    // Write data
    if (params_.write.value()->outputFileFormat.value() == "arome legacy binary"
      || params_.write.value()->outputFileFormat.value() == "arome legacy netcdf") {
      // Create AROME covariance fieldset
      atlas::FieldSet aromeCovData;

      // Compute covariance from correlation square-root and standard-deviation if it is missing
      covar_->computeCovariance(aromeCovData);

      // Define global vectors
      std::vector<double> vorCovGlb;
      std::vector<double> divuCovGlb;
      std::vector<double> tPsuCovGlb;
      std::vector<double> quCovGlb;

      for (const auto & var : centralVars()) {
        // Get number of levels
        const size_t nz = var.getLevels();

        // Get covariance view
        const auto covView = getView3D("cov", var, aromeCovData);

        // Create AROME covariance field
        createField3D("aromeCov", covar_->trans()->nw(), var, aromeCovData);

        // Get AROME covariance view
        auto aromeCovView = getView3D("aromeCov", var, aromeCovData);

        for (size_t jw = 0; jw < covar_->trans()->nw(); ++jw) {
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

      for (const auto & var : centralVars()) {
        // Get covariance field
        const auto aromeCovField = getField("aromeCov", var, aromeCovData);

        // Gather covariance vector
        if (var.name() == "air_upward_absolute_vorticity") {
          covar_->trans()->gatherCov(aromeCovField, vorCovGlb, true);
        }
        if (var.name() == "air_horizontal_divergence") {
          covar_->trans()->gatherCov(aromeCovField, divuCovGlb, true);
        }
        if (var.name() == "air_temperature_and_log_of_air_pressure_at_surface") {
          covar_->trans()->gatherCov(aromeCovField, tPsuCovGlb, true);
        }
        if (var.name() == "water_vapor_mixing_ratio_wrt_moist_air") {
          covar_->trans()->gatherCov(aromeCovField, quCovGlb, true);
        }
      }

      // Define attributes
      eckit::LocalConfiguration attributes;
      const atlas::StructuredGrid & outerGrid =
        covar_->trans()->geometryData().functionSpace().grid();
      const atlas::StructuredGrid & gpGrid = covar_->trans()->gpFspace().grid();
      const bool y_increasing = outerGrid.spec().getSubConfiguration("yspace").getDouble("end")
        > outerGrid.spec().getSubConfiguration("yspace").getDouble("start");
      attributes.set("clid", "ALADIN98");
      attributes.set("clcom", " Balanced statistcs for a LAM, after L. Berre 1998");
      attributes.set("iorig", 85);
      attributes.set("elon0", outerGrid.projection().spec().getDouble("longitude0"));
      attributes.set("elat0", outerGrid.projection().spec().getDouble("latitude0"));
      const auto corner1 = y_increasing ? outerGrid.lonlat(0, 0)
        : outerGrid.lonlat(0, outerGrid.ny()-1);
      attributes.set("elon1", corner1[0]);
      attributes.set("elat1", corner1[1]);
      const auto corner2 = y_increasing ? outerGrid.lonlat(outerGrid.nxmax()-1, outerGrid.ny()-1)
        :  outerGrid.lonlat(outerGrid.nxmax()-1, 0);
      attributes.set("elon2", corner2[0]);
      attributes.set("elat2", corner2[1]);
      attributes.set("ndgl", gpGrid.ny());
      attributes.set("ndlon", gpGrid.nxmax());
      attributes.set("ndgux", outerGrid.ny());
      attributes.set("ndlux", outerGrid.nxmax());
      attributes.set("nsmax", covar_->trans()->nwGlb()-1);
      attributes.set("nmsmax", covar_->trans()->M());
      attributes.set("nflev", centralVars()["air_upward_absolute_vorticity"].getLevels());

      if (covar_->comm().rank() == 0) {
        if (params_.write.value()->outputFileFormat.value() == "arome legacy binary") {
          // Write Fortran unformatted file (from ewgsacov.F90)
          const int nsmax = attributes.getDouble("nsmax");
          const int nflev = attributes.getDouble("nflev");
          bifourier_arome_legacy_write_covariance_f90(params_.write.value()->toConfiguration(),
            attributes, nsmax+1, nflev, vorCovGlb.data(), divuCovGlb.data(),
            tPsuCovGlb.data(), quCovGlb.data());
        } else if (params_.write.value()->outputFileFormat.value() == "arome legacy netcdf") {
          // NetCDF file path
          const std::string ncFilePath = params_.write.value()->outputFile.value();

          // NetCDF IDs
          int ncId, retval, nflevId, nflevp1Id, nsmaxp1Id, dNzId[3], dNzP1Id[3],
            vorCovId, divuCovId, tPsuCovId, quCovId;

          // Create NetCDF file
          if ((retval = nc_create(ncFilePath.c_str(), NC_64BIT_OFFSET | NC_CLOBBER, &ncId)))
            ERR(retval, ncFilePath);

          // Define attributes
          const std::string clid = attributes.getString("clid");
          if ((retval = nc_put_att_text(ncId, NC_GLOBAL, "ID", clid.size(), clid.c_str())))
            ERR(retval, "ID");
          const std::string clcom = attributes.getString("clcom");
          if ((retval = nc_put_att_text(ncId, NC_GLOBAL, "COMMENT", clcom.size(), clcom.c_str())))
            ERR(retval, "COMMENT");
          const int iorig = attributes.getInt("iorig");
          if ((retval = nc_put_att_int(ncId, NC_GLOBAL, "WMO_CENTRE", NC_INT, 1, &iorig)))
            ERR(retval, "WMO_CENTRE");
          const double elon0 = attributes.getDouble("elon0");
          if ((retval = nc_put_att_double(ncId, NC_GLOBAL, "ELON0", NC_DOUBLE, 1, &elon0)))
            ERR(retval, "ELON0");
          const double elat0 = attributes.getDouble("elat0");
          if ((retval = nc_put_att_double(ncId, NC_GLOBAL, "ELAT0", NC_DOUBLE, 1, &elat0)))
            ERR(retval, "ELAT0");
          const double elon1 = attributes.getDouble("elon1");
          if ((retval = nc_put_att_double(ncId, NC_GLOBAL, "ELON1", NC_DOUBLE, 1, &elon1)))
            ERR(retval, "ELON1");
          const double elat1 = attributes.getDouble("elat1");
          if ((retval = nc_put_att_double(ncId, NC_GLOBAL, "ELAT1", NC_DOUBLE, 1, &elat1)))
            ERR(retval, "ELAT1");
          const double elon2 = attributes.getDouble("elon2");
          if ((retval = nc_put_att_double(ncId, NC_GLOBAL, "ELON2", NC_DOUBLE, 1, &elon2)))
            ERR(retval, "ELON2");
          const double elat2 = attributes.getDouble("elat2");
          if ((retval = nc_put_att_double(ncId, NC_GLOBAL, "ELAT2", NC_DOUBLE, 1, &elat2)))
            ERR(retval, "ELAT2");
          const int ndgl = attributes.getDouble("ndgl");
          if ((retval = nc_put_att_int(ncId, NC_GLOBAL, "NDGL", NC_INT, 1, &ndgl)))
            ERR(retval, "NDGL");
          const int ndlon = attributes.getDouble("ndlon");
          if ((retval = nc_put_att_int(ncId, NC_GLOBAL, "NDLON", NC_INT, 1, &ndlon)))
            ERR(retval, "NDLON");
          const int ndgux = attributes.getDouble("ndgux");
          if ((retval = nc_put_att_int(ncId, NC_GLOBAL, "NDGUX", NC_INT, 1, &ndgux)))
            ERR(retval, "NDGUX");
          const int ndlux = attributes.getDouble("ndlux");
          if ((retval = nc_put_att_int(ncId, NC_GLOBAL, "NDLUX", NC_INT, 1, &ndlux)))
            ERR(retval, "NDLUX");
          const int nsmax = attributes.getDouble("nsmax");
          if ((retval = nc_put_att_int(ncId, NC_GLOBAL, "NSMAX", NC_INT, 1, &nsmax)))
            ERR(retval, "NSMAX");
          const int nmsmax = attributes.getDouble("nmsmax");
          if ((retval = nc_put_att_int(ncId, NC_GLOBAL, "NMSMAX", NC_INT, 1, &nmsmax)))
            ERR(retval, "NMSMAX");

          // Create dimensions
          const int nflev = attributes.getDouble("nflev");
          if ((retval = nc_def_dim(ncId, "NFLEV", nflev, &nflevId))) ERR(retval, "NFLEV");
          if ((retval = nc_def_dim(ncId, "NFLEVP1", nflev+1, &nflevp1Id))) ERR(retval, "NFLEVP1");
          if ((retval = nc_def_dim(ncId, "NSMAXP1", nsmax+1, &nsmaxp1Id)))
            ERR(retval, "NSMAXP1");

          // Dimensions arrays
          dNzId[0] = nsmaxp1Id;
          dNzId[1] = nflevId;
          dNzId[2] = nflevId;
          dNzP1Id[0] = nsmaxp1Id;
          dNzP1Id[1] = nflevp1Id;
          dNzP1Id[2] = nflevp1Id;

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
  const size_t jwGlb = jw + covar_->trans()->nwStart();

  // Constant coefficient
  const double zmovern = static_cast<double>(covar_->trans()->ellips().size())
    / static_cast<double>(covar_->trans()->nwGlb()-1);

  // Compute weight
  double zWeight;
  if (jwGlb != 0 && jwGlb != covar_->trans()->nwGlb()-1) {
    zWeight = 2.0*M_PI*static_cast<double>(jwGlb)*zmovern;
  } else if (jwGlb == 0) {
//    zWeight = M_PI*zmovern/4.0;
    zWeight = M_PI*zmovern/2.0;
  } else if (jwGlb == covar_->trans()->nwGlb()-1) {
    zWeight = M_PI*(static_cast<double>(covar_->trans()->nwGlb()-1)-0.25)*zmovern;
  }

  // REDNMC factor
  zWeight /= params_.rednmc.value()*params_.rednmc.value();

  oops::Log::trace() << classname() << "::aromeWeight starting" << std::endl;
  return zWeight;
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber

