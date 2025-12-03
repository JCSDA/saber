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

static SaberOuterBlockMaker<BifourierAromeBalance>
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
  ASSERT((params_.explicitPb.value() != boost::none) || params_.pbFromTrans.value()
    || params_.read.value());

  if ((params_.explicitPb.value() != boost::none) || params_.pbFromTrans.value()) {
    // Get change of variable parameters from configuration or from spectral transform
    const auto & explicitPb = params_.explicitPb.value();
    const size_t M = (explicitPb != boost::none) ? explicitPb->M.value() : trans_->M();
    const size_t N = (explicitPb != boost::none) ? explicitPb->N.value() : trans_->N();
    const double meanLat = (explicitPb != boost::none) ? explicitPb->meanLat.value()
      : trans_->meanLat();

    // Allocate fact1
    fact1_.resize(trans_->ns());

    // Compute change of variable factor
    const size_t nwGlb = std::max(M, N)+1;
    const double zromega = 0.7292115e-4;
    const double zcc = -2.0*zromega*std::sin(meanLat*M_PI/180.0);
    const double zly = 2.0*static_cast<double>(nwGlb)*trans_->dy();
    const double zfact1 = zcc*(zly/(2.0*M_PI))*(zly/(2.0*M_PI));
    for (size_t js = 0; js < trans_->ns(); ++js) {
      const double kstar = trans_->rkstar(trans_->jk(js), trans_->jl(js), M, N, nwGlb);
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
    size_t nial = 0;
    for (size_t jm = 0; jm < trans_->ellips().size(); ++jm) {
      nial += 4*(trans_->ellips()[jm]+1);
    }

    // Fact1 IAL vector
    std::vector<double> fact1IAL(nial);

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
        bifourier_arome_legacy_read_balance_f90(params_.read.value()->toConfiguration(),
          trans_->nwGlb(), nz_, sDivPbGlb.data(), sTpsPbGlb.data(), sTpsDivuGlb.data(),
          sQPbGlb.data(), sQDivuGlb.data(), sQTpsuGlb.data(), nial, fact1IAL.data());
      } else if (params_.read.value()->inputFileFormat.value() == "arome legacy netcdf") {
        // NetCDF file path
        const std::string ncFilePath = params_.read.value()->inputFile.value();

        // NetCDF IDs
        int ncId, retval, dimId, varId;
        size_t nzFromFile, nwGlbFromFile, nialFromFile;

        // Open NetCDF file
        if ((retval = nc_open(ncFilePath.c_str(), NC_NOWRITE, &ncId))) ERR(retval, ncFilePath);

        // Check dimensions
        if ((retval = nc_inq_dimid(ncId, "NFLEV", &dimId))) ERR(retval, "NFLEV");
        if ((retval = nc_inq_dimlen(ncId, dimId, &nzFromFile))) ERR(retval, "NFLEV");
        ASSERT(nzFromFile == nz_);
        if ((retval = nc_inq_dimid(ncId, "NSMAXP1", &dimId))) ERR(retval, "NSMAXP1");
        if ((retval = nc_inq_dimlen(ncId, dimId, &nwGlbFromFile))) ERR(retval, "NSMAXP1");
        ASSERT(nwGlbFromFile == trans_->nwGlb());
        if ((retval = nc_inq_dimid(ncId, "KSPEC2G", &dimId))) ERR(retval, "KSPEC2G");
        if ((retval = nc_inq_dimlen(ncId, dimId, &nialFromFile))) ERR(retval, "KSPEC2G");
        ASSERT(nialFromFile == nial);

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
    ASSERT(jIAL == nial);

    // Copy fact1
    for (size_t js = 0; js < trans_->ns(); ++js) {
      const size_t jk = trans_->jk(js);
      const size_t jl = trans_->jl(js);
      const size_t jq = trans_->jq(js);
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
  if (!((params_.explicitPb.value() != boost::none) || params_.pbFromTrans.value())) {
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

  if (params_.write.value() != boost::none) {
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
      size_t nial = 0;
      for (size_t jm = 0; jm < trans_->ellips().size(); ++jm) {
        nial += 4*(trans_->ellips()[jm]+1);
      }

      // Allocate fact1 IAL vector
      std::vector<double> fact1IAL(nial, 0.0);

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
      ASSERT(jIAL == nial);

      // Copy fact1
      for (size_t js = 0; js < trans_->ns(); ++js) {
        const size_t jk = trans_->jk(js);
        const size_t jl = trans_->jl(js);
        const size_t jq = trans_->jq(js);
        jIAL = IALIndexView(jk, jl, jq);
        fact1IAL[jIAL] = fact1_[js];
      }

      // Reduce fact1 IAL vector
      comm_.allReduceInPlace(fact1IAL.begin(), fact1IAL.end(), eckit::mpi::sum());

      if (comm_.rank() == 0) {
        // Get number of levels
        const size_t nz = balVars_["balanced_air_pressure"].getLevels();

        if (params_.write.value()->outputFileFormat.value() == "arome legacy binary") {
          // Write Fortran unformatted file (based on ewgsabal.F90)
          bifourier_arome_legacy_write_balance_f90(params_.write.value()->toConfiguration(),
            trans_->nwGlb(), nz_, sDivPbGlb.data(), sTpsPbGlb.data(), sTpsDivuGlb.data(),
            sQPbGlb.data(), sQDivuGlb.data(), sQTpsuGlb.data(), nial, fact1IAL.data());
        } else if (params_.write.value()->outputFileFormat.value() == "arome legacy netcdf") {
          // NetCDF file path
          const std::string ncFilePath = params_.write.value()->outputFile.value();

          // NetCDF IDs
          int ncId, retval, nzId, nzP1Id, nwGlbId, nialId, dNzNzId[3], dNzNzP1Id[3], dNzP1NzId[3],
            dIALId[1], sDivPbID, sTpsPbId, sTpsDivuId, sQPbId, sQDivuId, sQTpsuId, fact1Id;

          // Create NetCDF file
          if ((retval = nc_create(ncFilePath.c_str(), NC_64BIT_OFFSET | NC_CLOBBER, &ncId)))
            ERR(retval, ncFilePath);

          // Create dimensions
          if ((retval = nc_def_dim(ncId, "NFLEV", nz, &nzId))) ERR(retval, "NFLEV");
          if ((retval = nc_def_dim(ncId, "NFLEVP1", nz+1, &nzP1Id))) ERR(retval, "NFLEVP1");
          if ((retval = nc_def_dim(ncId, "NSMAXP1", trans_->nwGlb(), &nwGlbId)))
            ERR(retval, "NSMAXP1");
          if ((retval = nc_def_dim(ncId, "KSPEC2G", nial, &nialId))) ERR(retval, "KSPEC2G");

          // Dimensions arrays
          dNzNzId[0] = nwGlbId;
          dNzNzId[1] = nzId;
          dNzNzId[2] = nzId;
          dNzNzP1Id[0] = nwGlbId;
          dNzNzP1Id[1] = nzId;
          dNzNzP1Id[2] = nzP1Id;
          dNzP1NzId[0] = nwGlbId;
          dNzP1NzId[1] = nzP1Id;
          dNzP1NzId[2] = nzId;
          dIALId[0] = nialId;

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
  atlas::Field pbField = trans_->spFspace()->createField<double>(
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
  atlas::Field pbField = trans_->spFspace()->createField<double>(
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
  atlas::Field pbField = trans_->spFspace()->createField<double>(
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
  atlas::Field tField = trans_->spFspace()->createField<double>(
    atlas::option::name("air_temperature") | atlas::option::levels(nz_));
  atlas::Field psField = trans_->spFspace()->createField<double>(
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
  atlas::Field tPsField = trans_->spFspace()->createField<double>(
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
