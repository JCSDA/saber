/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/BifourierAromeBalance.h"

#include <netcdf.h>

#include <algorithm>

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

static OuterBlockMaker<BifourierAromeBalance>
  makerBifourierAromeBalance_("BifourierAromeBalance");

// -----------------------------------------------------------------------------

BifourierAromeBalance::BifourierAromeBalance(const oops::GeometryData & outerGeometryData,
                                   const oops::Variables & outerVars,
                                   const eckit::Configuration & covarConfig,
                                   const Parameters_ & params,
                                   const oops::FieldSet3D & xb,
                                   const oops::FieldSet3D & fg)
  : BifourierBalance(outerGeometryData, genericInnerVars(outerVars), covarConfig, params, xb, fg),
    params_(params),
    aromeInnerVars_(innerVars_)
{
  oops::Log::trace() << classname() << "::BifourierAromeBalance starting" << std::endl;

  // Check balanced air pressure source
  ASSERT((params_.explicitPb.value()) || params_.pbFromTrans.value()
    || params_.read.value());

  if ((params_.explicitPb.value()) || params_.pbFromTrans.value()) {
    // Get change of variable parameters from configuration or from spectral transform
    const auto & explicitPb = params_.explicitPb.value();
    const size_t M = (explicitPb) ? explicitPb->M.value() : trans_->M();
    const size_t N = (explicitPb) ? explicitPb->N.value() : trans_->N();
    const double meanLat = (explicitPb) ? explicitPb->meanLat.value() : trans_->meanLat();

    // Allocate fact1
    fact1_.resize(trans_->ns());

    // Compute change of variable factor
    const size_t nwGlb = std::max(M, N)+1;
    const double zromega = 0.7292115e-4;
    const double zcc = -2.0*zromega*std::sin(meanLat*M_PI/180.0);
    const double zly = 2.0*static_cast<double>(nwGlb)*trans_->dy();
    const double zfact1 = zcc*(zly/(2.0*M_PI))*(zly/(2.0*M_PI));
    for (size_t js = 0; js < trans_->ns(); ++js) {
      const double kstar = trans_->rkstar(trans_->k(js), trans_->l(js), M, N, nwGlb);
      if (kstar > 0.0) {
        fact1_[js] = zfact1/(kstar*kstar);
      } else {
        fact1_[js] = 0.0;
      }
    }
  }

  // Remove balanced pressure from inner variables
  aromeInnerVars_ -= aromeInnerVars_["balanced_air_pressure"];

  oops::Log::trace() << classname() << "::BifourierAromeBalance done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierAromeBalance::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  // Vorticity to balanced pressure
  vorToPb(fset);

  // Generic balance
  BifourierBalance::multiply(fset);

  // Remove balanced pressure
  removePb(fset);

  // Split TPs
  splitTPs(fset);

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierAromeBalance::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  // Split TPs, adjoint
  gatherTPs(fset);

  // Remove balanced pressure, adjoint
  removePbAD(fset);

  // Generic balance, adjoint
  BifourierBalance::multiplyAD(fset);

  // Vorticity to balanced pressure, adjoint
  vorToPbAD(fset);

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierAromeBalance::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;

  // Split TPs, left inverse
  gatherTPs(fset);

  // Remove balanced pressure, left inverse
  removePbLeftInverse(fset);

  // Generic balance, left inverse
  BifourierBalance::leftInverseMultiply(fset);

  // Vorticity to balanced pressure, left inverse
  vorToPbLeftInverse(fset);

  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierAromeBalance::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;

  // Allocate fact1
  std::vector<double> fact1FromFile(trans_->ns());

  // Read data
  if (params_.read.value()->inputFileFormat.value() == "arome legacy binary"
    || params_.read.value()->inputFileFormat.value() == "arome legacy netcdf") {
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

    // Define global vectors
    std::vector<double> sDivPbGlb;
    std::vector<double> sTpsPbGlb;
    std::vector<double> sTpsDivuGlb;
    std::vector<double> sQPbGlb;
    std::vector<double> sQDivuGlb;
    std::vector<double> sQTpsuGlb;

    // Define global IAL size
    size_t kspec2g = 0;
    for (size_t jm = 0; jm < trans_->ellips().size(); ++jm) {
      kspec2g += 4*(trans_->ellips()[jm]+1);
    }

    // Fact1 IAL vector
    std::vector<double> fact1IAL(kspec2g);

    // Define attributes
    eckit::LocalConfiguration attributes;
    attributes.set("nsmax", trans_->nwGlb()-1);
    attributes.set("nflev", nz_);
    attributes.set("kspec2g", kspec2g);

    if (comm_.rank() == 0) {
      // Allocate global vectors
      sDivPbGlb.resize(trans_->nwGlb()*nz_*nz_);
      sTpsPbGlb.resize(trans_->nwGlb()*nz_*(nz_+1));
      sTpsDivuGlb.resize(trans_->nwGlb()*nz_*(nz_+1));
      sQPbGlb.resize(trans_->nwGlb()*nz_*nz_);
      sQDivuGlb.resize(trans_->nwGlb()*nz_*nz_);
      sQTpsuGlb.resize(trans_->nwGlb()*(nz_+1)*nz_);

      if (params_.read.value()->inputFileFormat.value() == "arome legacy binary") {
        // Read Fortran unformatted file (based on readjbbal.F90)
        const int nsmax = attributes.getDouble("nsmax");
        const int nflev = attributes.getDouble("nflev");
        const int kspec2g = attributes.getDouble("kspec2g");
        bifourier_arome_legacy_read_balance_f90(params_.read.value()->toConfiguration(),
          attributes, nsmax+1, nflev, kspec2g, sDivPbGlb.data(), sTpsPbGlb.data(),
          sTpsDivuGlb.data(), sQPbGlb.data(), sQDivuGlb.data(), sQTpsuGlb.data(), fact1IAL.data());
      } else if (params_.read.value()->inputFileFormat.value() == "arome legacy netcdf") {
        // NetCDF file path
        const std::string ncFilePath = params_.read.value()->inputFile.value();

        // NetCDF IDs
        int ncId, retval, dimId, varId;
        size_t nflevFile, nsmaxp1File, kspec2gFile;

        // Open NetCDF file
        if ((retval = nc_open(ncFilePath.c_str(), NC_NOWRITE, &ncId))) ERR(retval, ncFilePath);

        // Check dimensions
        if ((retval = nc_inq_dimid(ncId, "NFLEV", &dimId))) ERR(retval, "NFLEV");
        if ((retval = nc_inq_dimlen(ncId, dimId, &nflevFile))) ERR(retval, "NFLEV");
        ASSERT(nflevFile == nz_);
        if ((retval = nc_inq_dimid(ncId, "NSMAXP1", &dimId))) ERR(retval, "NSMAXP1");
        if ((retval = nc_inq_dimlen(ncId, dimId, &nsmaxp1File))) ERR(retval, "NSMAXP1");
        ASSERT(nsmaxp1File == trans_->nwGlb());
        if ((retval = nc_inq_dimid(ncId, "KSPEC2G", &dimId))) ERR(retval, "KSPEC2G");
        if ((retval = nc_inq_dimlen(ncId, dimId, &kspec2gFile))) ERR(retval, "KSPEC2G");
        ASSERT(kspec2gFile == kspec2g);

        // Get variables
        if ((retval = nc_inq_varid(ncId, "SDIV_PB", &varId))) ERR(retval, "SDIV_PB");
        if ((retval = nc_get_var_double(ncId, varId, sDivPbGlb.data()))) ERR(retval, "SDIV_PB");
        if ((retval = nc_inq_varid(ncId, "STPS_PB", &varId))) ERR(retval, "STPS_PB");
        if ((retval = nc_get_var_double(ncId, varId, sTpsPbGlb.data()))) ERR(retval, "STPS_PB");
        if ((retval = nc_inq_varid(ncId, "STPS_DIVU", &varId))) ERR(retval, "STPS_DIVU");
        if ((retval = nc_get_var_double(ncId, varId, sTpsDivuGlb.data()))) ERR(retval, "STPS_DIVU");
        if ((retval = nc_inq_varid(ncId, "SQ_PB", &varId))) ERR(retval, "SQ_PB");
        if ((retval = nc_get_var_double(ncId, varId, sQPbGlb.data()))) ERR(retval, "SQ_PB");
        if ((retval = nc_inq_varid(ncId, "SQ_DIVU", &varId))) ERR(retval, "SQ_DIVU");
        if ((retval = nc_get_var_double(ncId, varId, sQDivuGlb.data()))) ERR(retval, "SQ_DIVU");
        if ((retval = nc_inq_varid(ncId, "SQ_TPSU", &varId))) ERR(retval, "SQ_TPSU");
        if ((retval = nc_get_var_double(ncId, varId, sQTpsuGlb.data()))) ERR(retval, "SQ_TPSU");
        if ((retval = nc_inq_varid(ncId, "FACT1", &varId))) ERR(retval, "FACT1");
        if ((retval = nc_get_var_double(ncId, varId, fact1IAL.data()))) ERR(retval, "FACT1");

        // Close file
        if ((retval = nc_close(ncId))) ERR(retval, ncFilePath);
      }
    }

    // Get fields
    auto sDivPbField = getField("reg", balVars_["air_horizontal_divergence"],
      balVars_["balanced_air_pressure"], data_);
    auto sTpsPbField = getField("reg",
      balVars_["air_temperature_and_log_of_air_pressure_at_surface"],
      balVars_["balanced_air_pressure"], data_);
    auto sTpsDivuField = getField("reg",
      balVars_["air_temperature_and_log_of_air_pressure_at_surface"],
      balVars_["air_horizontal_divergence"], data_);
    auto sQPbField = getField("reg", balVars_["water_vapor_mixing_ratio_wrt_moist_air"],
      balVars_["balanced_air_pressure"], data_);
    auto sQDivuField = getField("reg", balVars_["water_vapor_mixing_ratio_wrt_moist_air"],
      balVars_["air_horizontal_divergence"], data_);
    auto sQTpsuField = getField("reg", balVars_["water_vapor_mixing_ratio_wrt_moist_air"],
      balVars_["air_temperature_and_log_of_air_pressure_at_surface"], data_);

    // Scatter vectors
    trans_->scatterCov(sDivPbGlb, sDivPbField, true);
    trans_->scatterCov(sTpsPbGlb, sTpsPbField, true);
    trans_->scatterCov(sTpsDivuGlb, sTpsDivuField, true);
    trans_->scatterCov(sQPbGlb, sQPbField, true);
    trans_->scatterCov(sQDivuGlb, sQDivuField, true);
    trans_->scatterCov(sQTpsuGlb, sQTpsuField, true);

    // Broadcast fact1
    oops::Log::info() << "Info     : Broadcast fact1" << std::endl;
    comm_.broadcast(fact1IAL.begin(), fact1IAL.end(), 0);

    // Global IAL / spectral conversion
    atlas::Field IALIndexField("IALIndex", make_datatype<int>(),
      make_shape(trans_->ellips().size(), trans_->ellips()[0]+1, 4));
    auto IALIndexView = make_view<int, 3>(IALIndexField);
    IALIndexView.assign(-1);
    size_t jIAL = 0;
    for (size_t jk = 0; jk < trans_->ellips().size(); ++jk) {
      for (size_t jl = 0; jl <= trans_->ellips()[jk]; ++jl) {
        for (size_t jq = 0; jq < 4; ++jq) {
          IALIndexView(jk, jl, jq) = jIAL;
          ++jIAL;
        }
      }
    }
    ASSERT(jIAL == kspec2g);

    // Copy fact1
    for (size_t js = 0; js < trans_->ns(); ++js) {
      const size_t jk = trans_->k(js);
      const size_t jl = trans_->l(js);
      const size_t jq = trans_->q(js);
      jIAL = IALIndexView(jk, jl, jq);
      fact1FromFile[js] = fact1IAL[jIAL];
    }

    // Print norms
    print(oops::Log::test());
  } else {
    // Read generic balance
    BifourierBalance::read();

    // NetCDF file path
    const std::string ncFilePath = params_.read.value()->inputFile.value();

    // NetCDF IDs
    int ncId, retval, nsGlbId, varId;
    size_t nsGlbFromFile;

    // Define global vector
    std::vector<double> fact1Glb;

    if (comm_.rank() == 0) {
      // Open NetCDF file
      if ((retval = nc_open(ncFilePath.c_str(), NC_NOWRITE, &ncId))) ERR(retval, ncFilePath);

      // Check dimension
      if ((retval = nc_inq_dimid(ncId, "nsGlb", &nsGlbId))) ERR(retval, "nsGlb");
      if ((retval = nc_inq_dimlen(ncId, nsGlbId, &nsGlbFromFile))) ERR(retval, "nsGlb");
      ASSERT(nsGlbFromFile == trans_->nsGlb());

      // Get variable ID
      if ((retval = nc_inq_varid(ncId, "fact1", &varId))) ERR(retval, "fact1");

      // Read data
      std::vector<double> fact1GlbOrdered(trans_->nsGlb());
      if ((retval = nc_get_var_double(ncId, varId, fact1GlbOrdered.data()))) ERR(retval, "fact1");

      // Reorder data
      fact1Glb.resize(trans_->nsGlb());
      for (size_t jsGlb = 0; jsGlb < trans_->nsGlb(); ++jsGlb) {
        fact1Glb[jsGlb] = fact1GlbOrdered[trans_->sMapping()[jsGlb]];
      }

      // Close file
      if ((retval = nc_close(ncId))) ERR(retval, ncFilePath);
    }

    // Scatter vector
    comm_.scatterv(fact1Glb.cbegin(), fact1Glb.cend(), trans_->sCounts(), trans_->sDispls(),
      fact1FromFile.begin(), fact1FromFile.end(), 0);
  }

  // Copy fact1 from file if it has not been defined in the constructor
  if (!((params_.explicitPb.value()) || params_.pbFromTrans.value())) {
    // Allocate fact1
    fact1_.resize(trans_->ns());

    // Copy fact1
    fact1_ = fact1FromFile;
  }

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierAromeBalance::directCalibration(const oops::FieldSets & fsetEns) {
  oops::Log::trace() << classname() << "::directCalibration starting" << std::endl;

  // Copy ensemble
  auto fsetEnsCopy = fsetEns;

  for (size_t je = 0; je < fsetEnsCopy.size(); ++je) {
    // Split TPs, left inverse
    gatherTPs(fsetEnsCopy[je]);

    // Remove balanced pressure, left inverse
    removePbLeftInverse(fsetEnsCopy[je]);
  }

  // Generic balance
  BifourierBalance::directCalibration(fsetEnsCopy);

  oops::Log::trace() << classname() << "::directCalibration done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierAromeBalance::iterativeCalibrationUpdate(const oops::FieldSet3D & fset) {
  oops::Log::trace() << classname() << "::iterativeCalibrationUpdate starting" << std::endl;

  // Copy fieldset
  auto fsetCopy = fset;

  // Split TPs, left inverse
  gatherTPs(fsetCopy);

  // Remove balanced pressure, left inverse
  removePbLeftInverse(fsetCopy);

  // Generic balance
  BifourierBalance::iterativeCalibrationUpdate(fsetCopy);

  oops::Log::trace() << classname() << "::iterativeCalibrationUpdate done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierAromeBalance::write() const {
  oops::Log::trace() << classname() << "::write starting" << std::endl;

  if (params_.write.value()) {
    // Write data
    if (params_.write.value()->outputFileFormat.value() == "arome legacy binary"
      || params_.write.value()->outputFileFormat.value() == "arome legacy netcdf") {
      // Define global vectors
      std::vector<double> sDivPbGlb;
      std::vector<double> sTpsPbGlb;
      std::vector<double> sTpsDivuGlb;
      std::vector<double> sQPbGlb;
      std::vector<double> sQDivuGlb;
      std::vector<double> sQTpsuGlb;

      // Get fields
      const auto sDivPbField = getField("reg", balVars_["air_horizontal_divergence"],
        balVars_["balanced_air_pressure"], data_);
      const auto sTpsPbField = getField("reg",
        balVars_["air_temperature_and_log_of_air_pressure_at_surface"],
        balVars_["balanced_air_pressure"], data_);
      const auto sTpsDivuField = getField("reg",
        balVars_["air_temperature_and_log_of_air_pressure_at_surface"],
        balVars_["air_horizontal_divergence"], data_);
      const auto sQPbField = getField("reg", balVars_["water_vapor_mixing_ratio_wrt_moist_air"],
        balVars_["balanced_air_pressure"], data_);
      const auto sQDivuField = getField("reg", balVars_["water_vapor_mixing_ratio_wrt_moist_air"],
       balVars_["air_horizontal_divergence"], data_);
      const auto sQTpsuField = getField("reg", balVars_["water_vapor_mixing_ratio_wrt_moist_air"],
        balVars_["air_temperature_and_log_of_air_pressure_at_surface"], data_);

      // Gather vectors
      trans_->gatherCov(sDivPbField, sDivPbGlb, true);
      trans_->gatherCov(sTpsPbField, sTpsPbGlb, true);
      trans_->gatherCov(sTpsDivuField, sTpsDivuGlb, true);
      trans_->gatherCov(sQPbField, sQPbGlb, true);
      trans_->gatherCov(sQDivuField, sQDivuGlb, true);
      trans_->gatherCov(sQTpsuField, sQTpsuGlb, true);

      // Define global IAL size
      size_t kspec2g = 0;
      for (size_t jm = 0; jm < trans_->ellips().size(); ++jm) {
        kspec2g += 4*(trans_->ellips()[jm]+1);
      }

      // Allocate fact1 IAL vector
      std::vector<double> fact1IAL(kspec2g, 0.0);

      // Global IAL / spectral conversion
      atlas::Field IALIndexField("IALIndex", make_datatype<int>(),
        make_shape(trans_->ellips().size(), trans_->ellips()[0]+1, 4));
      auto IALIndexView = make_view<int, 3>(IALIndexField);
      IALIndexView.assign(-1);
      size_t jIAL = 0;
      for (size_t jk = 0; jk < trans_->ellips().size(); ++jk) {
        for (size_t jl = 0; jl <= trans_->ellips()[jk]; ++jl) {
          for (size_t jq = 0; jq < 4; ++jq) {
            IALIndexView(jk, jl, jq) = jIAL;
            ++jIAL;
          }
        }
      }
      ASSERT(jIAL == kspec2g);

      // Copy fact1
      for (size_t js = 0; js < trans_->ns(); ++js) {
        const size_t jk = trans_->k(js);
        const size_t jl = trans_->l(js);
        const size_t jq = trans_->q(js);
        jIAL = IALIndexView(jk, jl, jq);
        fact1IAL[jIAL] = fact1_[js];
      }

      // Reduce fact1 IAL vector
      comm_.allReduceInPlace(fact1IAL.begin(), fact1IAL.end(), eckit::mpi::sum());

      // Define attributes
      eckit::LocalConfiguration attributes;
      const atlas::StructuredGrid & outerGrid = trans_->geometryData().functionSpace().grid();
      const atlas::StructuredGrid & gpGrid = trans_->gpFspace().grid();
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
      attributes.set("nsmax", trans_->nwGlb()-1);
      attributes.set("nmsmax", trans_->M());
      attributes.set("nflev", nz_);
      attributes.set("kspec2g", kspec2g);

      if (comm_.rank() == 0) {
        if (params_.write.value()->outputFileFormat.value() == "arome legacy binary") {
          // Write Fortran unformatted file (based on ewgsabal.F90)
          const int nsmax = attributes.getDouble("nsmax");
          const int nflev = attributes.getDouble("nflev");
          const int kspec2g = attributes.getDouble("kspec2g");
          bifourier_arome_legacy_write_balance_f90(params_.write.value()->toConfiguration(),
            attributes, nsmax+1, nflev, kspec2g, sDivPbGlb.data(), sTpsPbGlb.data(),
            sTpsDivuGlb.data(), sQPbGlb.data(), sQDivuGlb.data(), sQTpsuGlb.data(),
            fact1IAL.data());
        } else if (params_.write.value()->outputFileFormat.value() == "arome legacy netcdf") {
          // NetCDF file path
          const std::string ncFilePath = params_.write.value()->outputFile.value();

          // NetCDF IDs
          int ncId, retval, nflevId, nflevp1Id, nsmaxp1Id, kspec2gId, dNzNzId[3], dNzNzP1Id[3],
            dNzP1NzId[3], dIALId[1], sDivPbID, sTpsPbId, sTpsDivuId, sQPbId, sQDivuId, sQTpsuId,
            fact1Id;

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
          const int kspec2g = attributes.getDouble("kspec2g");
          if ((retval = nc_def_dim(ncId, "KSPEC2G", kspec2g, &kspec2gId))) ERR(retval, "KSPEC2G");

          // Dimensions arrays
          dNzNzId[0] = nsmaxp1Id;
          dNzNzId[1] = nflevId;
          dNzNzId[2] = nflevId;
          dNzNzP1Id[0] = nsmaxp1Id;
          dNzNzP1Id[1] = nflevId;
          dNzNzP1Id[2] = nflevp1Id;
          dNzP1NzId[0] = nsmaxp1Id;
          dNzP1NzId[1] = nflevp1Id;
          dNzP1NzId[2] = nflevId;
          dIALId[0] = kspec2gId;

          // Create variables
          if ((retval = nc_def_var(ncId, "SDIV_PB", NC_DOUBLE, 3, dNzNzId, &sDivPbID)))
            ERR(retval, "SDIV_PB");
          if ((retval = nc_def_var(ncId, "STPS_PB", NC_DOUBLE, 3, dNzNzP1Id, &sTpsPbId)))
            ERR(retval, "STPS_PB");
          if ((retval = nc_def_var(ncId, "STPS_DIVU", NC_DOUBLE, 3, dNzNzP1Id, &sTpsDivuId)))
            ERR(retval, "STPS_DIVU");
          if ((retval = nc_def_var(ncId, "SQ_PB", NC_DOUBLE, 3, dNzNzId, &sQPbId)))
            ERR(retval, "SQ_PB");
          if ((retval = nc_def_var(ncId, "SQ_DIVU", NC_DOUBLE, 3, dNzNzId, &sQDivuId)))
            ERR(retval, "SQ_DIVU");
          if ((retval = nc_def_var(ncId, "SQ_TPSU", NC_DOUBLE, 3, dNzP1NzId, &sQTpsuId)))
            ERR(retval, "SQ_TPSU");
          if ((retval = nc_def_var(ncId, "FACT1", NC_DOUBLE, 1, dIALId, &fact1Id)))
            ERR(retval, "FACT1");

          // End definition mode
          if ((retval = nc_enddef(ncId))) ERR(retval, ncFilePath);

          // Write data
          if ((retval = nc_put_var_double(ncId, sDivPbID, sDivPbGlb.data())))
            ERR(retval, "SDIV_PB");
          if ((retval = nc_put_var_double(ncId, sTpsPbId, sTpsPbGlb.data())))
            ERR(retval, "STPS_PB");
          if ((retval = nc_put_var_double(ncId, sTpsDivuId, sTpsDivuGlb.data())))
            ERR(retval, "STPS_DIVU");
          if ((retval = nc_put_var_double(ncId, sQPbId, sQPbGlb.data())))
            ERR(retval, "SQ_PB");
          if ((retval = nc_put_var_double(ncId, sQDivuId, sQDivuGlb.data())))
            ERR(retval, "SQ_DIVU");
          if ((retval = nc_put_var_double(ncId, sQTpsuId, sQTpsuGlb.data())))
            ERR(retval, "SQ_TPSU");
          if ((retval = nc_put_var_double(ncId, fact1Id, fact1IAL.data())))
            ERR(retval, "FACT1");

          // Close file
          if ((retval = nc_close(ncId))) ERR(retval, ncFilePath);
        }
      }
    } else {
      // Generic balance
      BifourierBalance::write();

      // Allocate global vector
      std::vector<double> fact1Glb;
      if (comm_.rank() == 0) {
        fact1Glb.resize(trans_->nsGlb());
      }

      // Gather data
      comm_.gatherv(fact1_.cbegin(), fact1_.cend(), fact1Glb.begin(), fact1Glb.end(),
        trans_->sCounts(), trans_->sDispls(), 0);

      // NetCDF IDs
      int retval, ncId, nsGlbId, d1DId[1], varId;

      // NetCDF file path
      const std::string ncFilePath = params_.write.value()->outputFile.value();

      if (comm_.rank() == 0) {
        // Open NetCDF file
        if ((retval = nc_open(ncFilePath.c_str(), NC_64BIT_OFFSET | NC_WRITE, &ncId)))
          ERR(retval, ncFilePath);

        // Return to definition mode
        if ((retval = nc_redef(ncId))) ERR(retval, ncFilePath);

        // Create dimension
        if ((retval = nc_def_dim(ncId, "nsGlb", trans_->nsGlb(), &nsGlbId))) ERR(retval, "nsGlb");

        // Dimensions array
        d1DId[0] = nsGlbId;

        // Define variable
        if ((retval = nc_def_var(ncId, "fact1", NC_DOUBLE, 1, d1DId, &varId)))
          ERR(retval, "fact1");

        // End definition mode
        if ((retval = nc_enddef(ncId))) ERR(retval, ncFilePath);

        // Reorder data
        std::vector<double> fact1GlbOrdered(trans_->nsGlb());
        for (size_t jsGlb = 0; jsGlb < trans_->nsGlb(); ++jsGlb) {
          fact1GlbOrdered[trans_->sMapping()[jsGlb]] = fact1Glb[jsGlb];
        }

        // Write data
        if ((retval = nc_put_var_double(ncId, varId, fact1GlbOrdered.data())))
          ERR(retval, "fact1");

        // Close file
        if ((retval = nc_close(ncId))) ERR(retval, ncFilePath);
      }
    }
  }

  oops::Log::trace() << classname() << "::write done" << std::endl;
}

// -----------------------------------------------------------------------------

oops::Variables BifourierAromeBalance::genericInnerVars(const oops::Variables & outerVars) {
  oops::Log::trace() << classname() << "::genericInnerVars starting" << std::endl;

  // Get number of levels
  nz_ = outerVars["air_temperature"].getLevels();

  // Add TPs to inner variables and remove T and Ps
  oops::Variables vars(outerVars);
  vars.push_back("air_temperature_and_log_of_air_pressure_at_surface");
  vars["air_temperature_and_log_of_air_pressure_at_surface"].setLevels(nz_+1);
  vars -= vars["air_temperature"];
  vars -= vars["log_of_air_pressure_at_surface"];

  // Add balanced pressure
  vars.push_back("balanced_air_pressure");
  vars["balanced_air_pressure"].setLevels(nz_);

  oops::Log::trace() << classname() << "::genericInnerVars done" << std::endl;
  return vars;
}

// -----------------------------------------------------------------------------

void BifourierAromeBalance::vorToPb(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::vorToPb starting" << std::endl;

  // Get inner field
  const auto vorField = fset["air_upward_absolute_vorticity"];

  // Create outer field
  atlas::Field pbField = trans_->spFspace().createField<double>(
    atlas::option::name("balanced_air_pressure") | atlas::option::levels(nz_));

  // Get fields views
  const auto vorView = make_view<double, 2>(vorField);
  auto pbView = make_view<double, 2>(pbField);

  // Apply change of variable
  for (size_t js = 0; js < trans_->ns(); ++js) {
    for (size_t jz = 0; jz < nz_; ++jz) {
      pbView(js, jz) = vorView(js, jz)*fact1_[js];
    }
  }

  // Add outer field
  fset.add(pbField);

  oops::Log::trace() << classname() << "::vorToPb done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierAromeBalance::vorToPbAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::vorToPbAD starting" << std::endl;

  // Get fields
  const auto pbField = fset["balanced_air_pressure"];
  auto vorField = fset["air_upward_absolute_vorticity"];

  // Get fields views
  const auto pbView = make_view<double, 2>(pbField);
  auto vorView = make_view<double, 2>(vorField);

  // Apply change of variable, adjoint
  for (size_t js = 0; js < trans_->ns(); ++js) {
    for (size_t jz = 0; jz < nz_; ++jz) {
      vorView(js, jz) += pbView(js, jz)*fact1_[js];
    }
  }

  // Remove outer field
  util::removeFieldsFromFieldSet(fset.fieldSet(), {"balanced_air_pressure"});

  oops::Log::trace() << classname() << "::vorToPbAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierAromeBalance::vorToPbLeftInverse(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::vorToPbLeftInverse starting" << std::endl;

  // Get fields
  const auto pbField = fset["balanced_air_pressure"];
  auto vorField = fset["air_upward_absolute_vorticity"];

  // Get fields views
  const auto pbView = make_view<double, 2>(pbField);
  auto vorView = make_view<double, 2>(vorField);

  // Apply change of variable, inverse
  for (size_t js = 0; js < trans_->ns(); ++js) {
    for (size_t jz = 0; jz < nz_; ++jz) {
      if (std::abs(fact1_[js]) > 0.0) {
        vorView(js, jz) = pbView(js, jz)/fact1_[js];
      }
    }
  }

  // Remove outer field
  util::removeFieldsFromFieldSet(fset.fieldSet(), {"balanced_air_pressure"});

  oops::Log::trace() << classname() << "::vorToPbLeftInverse done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierAromeBalance::removePb(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::removePb starting" << std::endl;

  util::removeFieldsFromFieldSet(fset.fieldSet(), {"balanced_air_pressure"});

  oops::Log::trace() << classname() << "::removePb done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierAromeBalance::removePbAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::removePbAD starting" << std::endl;

  // Get outer field
  const auto vorField = fset["air_upward_absolute_vorticity"];

  // Create inner field
  atlas::Field pbField = trans_->spFspace().createField<double>(
    atlas::option::name("balanced_air_pressure") | atlas::option::levels(nz_));

  // Get inner field view
  auto pbView = make_view<double, 2>(pbField);

  // Set inner field to zero
  pbView.assign(0.0);

  // Add outer field
  fset.add(pbField);

  oops::Log::trace() << classname() << "::removePbAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierAromeBalance::removePbLeftInverse(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::removePbLeftInverse starting" << std::endl;

  // Get inner field
  const auto vorField = fset["air_upward_absolute_vorticity"];

  // Create outer field
  atlas::Field pbField = trans_->spFspace().createField<double>(
    atlas::option::name("balanced_air_pressure") | atlas::option::levels(nz_));

  // Get fields views
  const auto vorView = make_view<double, 2>(vorField);
  auto pbView = make_view<double, 2>(pbField);

  // Apply change of variable
  for (size_t js = 0; js < trans_->ns(); ++js) {
    for (size_t jz = 0; jz < nz_; ++jz) {
      pbView(js, jz) = vorView(js, jz)*fact1_[js];
    }
  }

  // Add outer field
  fset.add(pbField);

  oops::Log::trace() << classname() << "::removePbLeftInverse done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierAromeBalance::splitTPs(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::splitTPs starting" << std::endl;

  // Get inner field
  const auto tPsField = fset["air_temperature_and_log_of_air_pressure_at_surface"];

  // Create outer fields
  atlas::Field tField = trans_->spFspace().createField<double>(
    atlas::option::name("air_temperature") | atlas::option::levels(nz_));
  atlas::Field psField = trans_->spFspace().createField<double>(
    atlas::option::name("log_of_air_pressure_at_surface") | atlas::option::levels(1));

  // Get fields views
  const auto tPsView = make_view<double, 2>(tPsField);
  auto tView = make_view<double, 2>(tField);
  auto psView = make_view<double, 2>(psField);

  // Copy data
  for (size_t js = 0; js < trans_->ns(); ++js) {
    for (size_t jz = 0; jz < nz_; ++jz) {
      tView(js, jz) = tPsView(js, jz);
    }
    psView(js, 0) = tPsView(js, nz_);
  }

  // Remove inner field
  util::removeFieldsFromFieldSet(fset.fieldSet(),
    {"air_temperature_and_log_of_air_pressure_at_surface"});

  // Add outer fields
  fset.add(tField);
  fset.add(psField);

  oops::Log::trace() << classname() << "::splitTPs done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierAromeBalance::gatherTPs(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::gatherTPs starting" << std::endl;

  // Get outer fields
  const auto tField = fset["air_temperature"];
  const auto psField = fset["log_of_air_pressure_at_surface"];

  // Create inner field
  atlas::Field tPsField = trans_->spFspace().createField<double>(
    atlas::option::name("air_temperature_and_log_of_air_pressure_at_surface") |
    atlas::option::levels(nz_+1));

  // Get fields views
  const auto tView = make_view<double, 2>(tField);
  const auto psView = make_view<double, 2>(psField);
  auto tPsView = make_view<double, 2>(tPsField);

  // Copy data
  for (size_t js = 0; js < trans_->ns(); ++js) {
    for (size_t jz = 0; jz < nz_; ++jz) {
      tPsView(js, jz) = tView(js, jz);
    }
    tPsView(js, nz_) = psView(js, 0);
  }

  // Remove outer fields
  util::removeFieldsFromFieldSet(fset.fieldSet(), {"air_temperature",
    "log_of_air_pressure_at_surface"});

  // Add inner field
  fset.add(tPsField);

  oops::Log::trace() << classname() << "::gatherTPs done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
