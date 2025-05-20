/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/BifourierAromeLegacy.h"

#include <Eigen/Dense>
#include <netcdf.h>

#include <algorithm>
#include <string>
#include <vector>

#include "oops/util/FieldSetHelpers.h"

#include "saber/bifourier/bifourier_arome_legacy.h"
#include "saber/bifourier/BifourierBalance.h"
#include "saber/bifourier/BifourierCovariance.h"
#include "saber/bifourier/BifourierUtilities.h"
#include "saber/bifourier/BifourierVorToPb.h"

#define ERR(e, msg) {std::string s(nc_strerror(e)); \
  throw eckit::Exception(s + " : " + msg, Here());}

using atlas::array::make_datatype;
using atlas::array::make_shape;
using atlas::array::make_view;

namespace saber {
namespace bifourier {
namespace arome_legacy {

// -----------------------------------------------------------------------------

void readVorToPb(const eckit::mpi::Comm & comm,
                 const BifourierVorToPbReadParameters & params,
                 const BifourierTransform & trans,
                 std::vector<double> & fact1) {
  oops::Log::trace() << "saber::bifourier::BifourierAromeLegacy::readVorToPb starting" << std::endl;

  // Define global IAL size
  size_t nial = 0;
  for (size_t jm = 0; jm < trans.ellips().size(); ++jm) {
    nial += 4*(trans.ellips()[jm]+1);
  }

  // Fact1 global vector
  std::vector<double> fact1Glb(nial);

  if (comm.rank() == 0) {
    if (params.inputFileFormat.value() == "arome legacy binary") {
      // Read Fortran unformatted file (based on readjbbal.F90)
      bifourier_arome_legacy_vortopb_f90(params.toConfiguration(), nial, fact1Glb.data());
    } else if (params.inputFileFormat.value() == "arome legacy netcdf") {
      // NetCDF file path
      const std::string ncFilePath = *params.inputFile.value();

      // NetCDF IDs
      int ncid, retval, dimid, varid;
      size_t nialFromFile;

      // Open NetCDF file
      if ((retval = nc_open(ncFilePath.c_str(), NC_NOWRITE, &ncid))) ERR(retval, ncFilePath);

      // Check dimensions
      if ((retval = nc_inq_dimid(ncid, "KSPEC2G", &dimid))) ERR(retval, "KSPEC2G");
      if ((retval = nc_inq_dimlen(ncid, dimid, &nialFromFile))) ERR(retval, "KSPEC2G");
      ASSERT(nialFromFile == nial);

      // Get variables
      if ((retval = nc_inq_varid(ncid, "FACT1", &varid))) ERR(retval, "FACT1");
      if ((retval = nc_get_var_double(ncid, varid, fact1Glb.data()))) ERR(retval, "FACT1");

      // Close file
      if ((retval = nc_close(ncid))) ERR(retval, ncFilePath);
    }

    // Print norms
    oops::Log::info() << "Info     : Fact1 norm: " << std::endl;
    const double fact1Norm = std::inner_product(fact1Glb.begin(), fact1Glb.end(), fact1Glb.begin(),
      0.0);
    oops::Log::test() << "- fact1: " << fact1Norm << std::endl;
  }

  // Broadcast fact1
  oops::Log::info() << "Info     : Broadcast fact1" << std::endl;
  comm.broadcast(fact1Glb.begin(), fact1Glb.end(), 0);

  // Global IAL / spectral conversion
  atlas::Field IALIndexField("IALIndex", make_datatype<int>(),
    make_shape(trans.ellips().size(), trans.ellips()[0]+1, 4));
  auto IALIndexView = make_view<int, 3>(IALIndexField);
  IALIndexView.assign(-1);
  size_t jIAL = 0;
  for (size_t jk = 0; jk < trans.ellips().size(); ++jk) {
    for (size_t jl = 0; jl <= trans.ellips()[jk]; ++jl) {
      for (size_t jq = 0; jq < 4; ++jq) {
        IALIndexView(jk, jl, jq) = jIAL;
        ++jIAL;
      }
    }
  }
  ASSERT(jIAL == nial);

  // Copy fact1
  for (size_t js = 0; js < trans.ns(); ++js) {
    const size_t jk = trans.jk(js);
    const size_t jl = trans.jl(js);
    const size_t jq = trans.jq(js);
    jIAL = IALIndexView(jk, jl, jq);
    fact1[js] = fact1Glb[jIAL];
  }

  oops::Log::trace() << "saber::bifourier::BifourierAromeLegacy::readVorToPb done" << std::endl;
}

// -----------------------------------------------------------------------------

void readBalance(const eckit::mpi::Comm & comm,
                 const oops::Variables & balVars,
                 const BifourierBalanceReadParameters & params,
                 const BifourierTransform & trans,
                 atlas::FieldSet & balFset) {
  oops::Log::trace() << "saber::bifourier::BifourierAromeLegacy::readBalance starting" << std::endl;

  // Get number of levels
  const size_t nflev = balVars[0].getLevels();

  // Allocate local vectors
  std::vector<double> sDivPb(trans.nw()*nflev*nflev);
  std::vector<double> sTpsPb(trans.nw()*nflev*(nflev+1));
  std::vector<double> sTpsDivu(trans.nw()*nflev*(nflev+1));
  std::vector<double> sQPb(trans.nw()*nflev*nflev);
  std::vector<double> sQDivu(trans.nw()*nflev*nflev);
  std::vector<double> sQTpsu(trans.nw()*(nflev+1)*nflev);

  // Tag root
  const size_t tagRoot = 123;

  if (comm.rank() == 0) {
    // Allocate global vectors
    std::vector<double> sDivPbGlb(trans.nwGlb()*nflev*nflev);
    std::vector<double> sTpsPbGlb(trans.nwGlb()*nflev*(nflev+1));
    std::vector<double> sTpsDivuGlb(trans.nwGlb()*nflev*(nflev+1));
    std::vector<double> sQPbGlb(trans.nwGlb()*nflev*nflev);
    std::vector<double> sQDivuGlb(trans.nwGlb()*nflev*nflev);
    std::vector<double> sQTpsuGlb(trans.nwGlb()*(nflev+1)*nflev);

    if (params.inputFileFormat.value() == "arome legacy binary") {
      // Read Fortran unformatted file (based on readjbbal.F90)
      bifourier_arome_legacy_balance_f90(params.toConfiguration(), trans.nwGlb(), nflev,
        sDivPbGlb.data(), sTpsPbGlb.data(), sTpsDivuGlb.data(), sQPbGlb.data(), sQDivuGlb.data(),
        sQTpsuGlb.data());
    } else if (params.inputFileFormat.value() == "arome legacy netcdf") {
      // NetCDF file path
      std::string ncFilePath = params.inputFile.value();

      // NetCDF IDs
      int ncid, retval, dimid, varid;
      size_t nflevFromFile, nwFromFile;

      // Open NetCDF file
      if ((retval = nc_open(ncFilePath.c_str(), NC_NOWRITE, &ncid))) ERR(retval, ncFilePath);

      // Check dimensions
      if ((retval = nc_inq_dimid(ncid, "NFLEV", &dimid))) ERR(retval, "NFLEV");
      if ((retval = nc_inq_dimlen(ncid, dimid, &nflevFromFile))) ERR(retval, "NFLEV");
      ASSERT(nflevFromFile == nflev);
      if ((retval = nc_inq_dimid(ncid, "NSMAXP1", &dimid))) ERR(retval, "NSMAXP1");
      if ((retval = nc_inq_dimlen(ncid, dimid, &nwFromFile))) ERR(retval, "NSMAXP1");
      ASSERT(nwFromFile == trans.nwGlb());

      // Get variables
      if ((retval = nc_inq_varid(ncid, "SDIV_PB", &varid))) ERR(retval, "SDIV_PB");
      if ((retval = nc_get_var_double(ncid, varid, sDivPbGlb.data()))) ERR(retval, "SDIV_PB");
      if ((retval = nc_inq_varid(ncid, "STPS_PB", &varid))) ERR(retval, "STPS_PB");
      if ((retval = nc_get_var_double(ncid, varid, sTpsPbGlb.data()))) ERR(retval, "STPS_PB");
      if ((retval = nc_inq_varid(ncid, "STPS_DIVU", &varid))) ERR(retval, "STPS_DIVU");
      if ((retval = nc_get_var_double(ncid, varid, sTpsDivuGlb.data()))) ERR(retval, "STPS_DIVU");
      if ((retval = nc_inq_varid(ncid, "SQ_PB", &varid))) ERR(retval, "SQ_PB");
      if ((retval = nc_get_var_double(ncid, varid, sQPbGlb.data()))) ERR(retval, "SQ_PB");
      if ((retval = nc_inq_varid(ncid, "SQ_DIVU", &varid))) ERR(retval, "SQ_DIVU");
      if ((retval = nc_get_var_double(ncid, varid, sQDivuGlb.data()))) ERR(retval, "SQ_DIVU");
      if ((retval = nc_inq_varid(ncid, "SQ_TPSU", &varid))) ERR(retval, "SQ_TPSU");
      if ((retval = nc_get_var_double(ncid, varid, sQTpsuGlb.data()))) ERR(retval, "SQ_TPSU");

      // Close file
      if ((retval = nc_close(ncid))) ERR(retval, ncFilePath);
    }

    // Print norms
    oops::Log::info() << "Info     : Balance operator norms: " << std::endl;
    const double sDivPbNorm = std::inner_product(sDivPbGlb.begin(), sDivPbGlb.end(),
      sDivPbGlb.begin(), 0.0);
    oops::Log::test() << "- sDivPb: " << sDivPbNorm << std::endl;
    const double sTpsPbNorm = std::inner_product(sTpsPbGlb.begin(), sTpsPbGlb.end(),
      sTpsPbGlb.begin(), 0.0);
    oops::Log::test() << "- sTpsPb: " << sTpsPbNorm << std::endl;
    const double sTpsDivuNorm = std::inner_product(sTpsDivuGlb.begin(), sTpsDivuGlb.end(),
      sTpsDivuGlb.begin(), 0.0);
    oops::Log::test() << "- sTpsDivu: " << sTpsDivuNorm << std::endl;
    const double sQPbNorm = std::inner_product(sQPbGlb.begin(), sQPbGlb.end(), sQPbGlb.begin(),
      0.0);
    oops::Log::test() << "- sQPb: " << sQPbNorm << std::endl;
    const double sQDivuNorm = std::inner_product(sQDivuGlb.begin(), sQDivuGlb.end(),
      sQDivuGlb.begin(), 0.0);
    oops::Log::test() << "- sQDivu: " << sQDivuNorm << std::endl;
    const double sQTpsuNorm = std::inner_product(sQTpsuGlb.begin(), sQTpsuGlb.end(),
      sQTpsuGlb.begin(), 0.0);
    oops::Log::test() << "- sQTpsu: " << sQTpsuNorm << std::endl;

    // Copy data
    for (size_t jw = 0; jw < trans.nw(); ++jw) {
      const size_t jwGlb = jw + trans.nwStart();
      for (size_t jz2 = 0; jz2 < nflev; ++jz2) {
        for (size_t jz1 = 0; jz1 < nflev; ++jz1) {
          const size_t index = jw*nflev*nflev + jz1*nflev + jz2;
          const size_t indexGlb = jwGlb*nflev*nflev + jz1*nflev + jz2;
          sDivPb[index] = sDivPbGlb[indexGlb];
          sQPb[index] = sQPbGlb[indexGlb];
          sQDivu[index] = sQDivuGlb[indexGlb];
        }
      }
      for (size_t jz2 = 0; jz2 < nflev+1; ++jz2) {
        for (size_t jz1 = 0; jz1 < nflev; ++jz1) {
          const size_t index = jw*nflev*(nflev+1) + jz1*(nflev+1) + jz2;
          const size_t indexGlb = jwGlb*nflev*(nflev+1) + jz1*(nflev+1) + jz2;
          sTpsPb[index] = sTpsPbGlb[indexGlb];
          sTpsDivu[index] = sTpsDivuGlb[indexGlb];
        }
      }
      for (size_t jz2 = 0; jz2 < nflev; ++jz2) {
        for (size_t jz1 = 0; jz1 < nflev+1; ++jz1) {
          const size_t index = jw*(nflev+1)*nflev + jz1*nflev + jz2;
          const size_t indexGlb = jwGlb*(nflev+1)*nflev + jz1*nflev + jz2;
          sQTpsu[index] = sQTpsuGlb[indexGlb];
        }
      }
    }

    // Send data
    int tag = tagRoot;
    for (size_t jt = 1; jt < comm.size(); ++jt) {
      // Create vectors
      std::vector<double> sDivPbSend(trans.nwPerTask()[jt]*nflev*nflev);
      std::vector<double> sTpsPbSend(trans.nwPerTask()[jt]*nflev*(nflev+1));
      std::vector<double> sTpsDivuSend(trans.nwPerTask()[jt]*nflev*(nflev+1));
      std::vector<double> sQPbSend(trans.nwPerTask()[jt]*nflev*nflev);
      std::vector<double> sQDivuSend(trans.nwPerTask()[jt]*nflev*nflev);
      std::vector<double> sQTpsuSend(trans.nwPerTask()[jt]*(nflev+1)*nflev);

      // Fill vectors
      for (size_t jwSend = 0; jwSend < trans.nwPerTask()[jt]; ++jwSend) {
        const size_t jwGlb = jwSend + trans.nwStartPerTask()[jt];
        for (size_t jz2 = 0; jz2 < nflev; ++jz2) {
          for (size_t jz1 = 0; jz1 < nflev; ++jz1) {
            const size_t indexSend = jwSend*nflev*nflev + jz1*nflev + jz2;
            const size_t indexGlb = jwGlb*nflev*nflev + jz1*nflev + jz2;
            sDivPbSend[indexSend] = sDivPbGlb[indexGlb];
            sQPbSend[indexSend] = sQPbGlb[indexGlb];
            sQDivuSend[indexSend] = sQDivuGlb[indexGlb];
          }
        }
        for (size_t jz2 = 0; jz2 < nflev+1; ++jz2) {
          for (size_t jz1 = 0; jz1 < nflev; ++jz1) {
            const size_t indexSend = jwSend*nflev*(nflev+1) + jz1*(nflev+1) + jz2;
            const size_t indexGlb = jwGlb*nflev*(nflev+1) + jz1*(nflev+1) + jz2;
            sTpsPbSend[indexSend] = sTpsPbGlb[indexGlb];
            sTpsDivuSend[indexSend] = sTpsDivuGlb[indexGlb];
          }
        }
        for (size_t jz2 = 0; jz2 < nflev; ++jz2) {
          for (size_t jz1 = 0; jz1 < nflev+1; ++jz1) {
            const size_t indexSend = jwSend*(nflev+1)*nflev + jz1*nflev + jz2;
            const size_t indexGlb = jwGlb*(nflev+1)*nflev + jz1*nflev + jz2;
            sQTpsuSend[indexSend] = sQTpsuGlb[indexGlb];
          }
        }
      }

      // Send vectors
      comm.send(sDivPbSend.data(), trans.nwPerTask()[jt]*nflev*nflev, jt, tag+0);
      comm.send(sTpsPbSend.data(), trans.nwPerTask()[jt]*nflev*(nflev+1), jt, tag+1);
      comm.send(sTpsDivuSend.data(), trans.nwPerTask()[jt]*nflev*(nflev+1), jt, tag+2);
      comm.send(sQPbSend.data(), trans.nwPerTask()[jt]*nflev*nflev, jt, tag+3);
      comm.send(sQDivuSend.data(), trans.nwPerTask()[jt]*nflev*nflev, jt, tag+4);
      comm.send(sQTpsuSend.data(), trans.nwPerTask()[jt]*(nflev+1)*nflev, jt, tag+5);

      // Update tag
      tag += 6;
    }
  }

  // Receive vectors
  if (comm.rank() > 0) {
    const int tagBase = tagRoot+6*(comm.rank()-1);
    comm.receive(sDivPb.data(), trans.nw()*nflev*nflev, 0, tagBase+0);
    comm.receive(sTpsPb.data(), trans.nw()*nflev*(nflev+1), 0, tagBase+1);
    comm.receive(sTpsDivu.data(), trans.nw()*nflev*(nflev+1), 0, tagBase+2);
    comm.receive(sQPb.data(), trans.nw()*nflev*nflev, 0, tagBase+3);
    comm.receive(sQDivu.data(), trans.nw()*nflev*nflev, 0, tagBase+4);
    comm.receive(sQTpsu.data(), trans.nw()*(nflev+1)*nflev, 0, tagBase+5);
  }

  // MPI barrier
  comm.barrier();

  // Get fields
  atlas::Field sDivPbField = balFset[
    "reg_air_horizontal_divergence_from_balanced_air_pressure"];
  atlas::Field sTpsPbField = balFset[
    "reg_air_temperature_and_log_of_air_pressure_at_surface_from_balanced_air_pressure"];
  atlas::Field sTpsDivuField = balFset[
    "reg_air_temperature_and_log_of_air_pressure_at_surface_from_air_horizontal_divergence"];
  atlas::Field sQPbField = balFset[
    "reg_water_vapor_mixing_ratio_wrt_moist_air_from_balanced_air_pressure"];
  atlas::Field sQDivuField = balFset[
    "reg_water_vapor_mixing_ratio_wrt_moist_air_from_air_horizontal_divergence"];
  atlas::Field sQTpsuField = balFset[
    "reg_water_vapor_mixing_ratio_wrt_moist_air_from_air_temperature_and_"
    "log_of_air_pressure_at_surface"];

  // Get fields views
  auto sDivPbView = make_view<double, 3>(sDivPbField);
  auto sTpsPbView = make_view<double, 3>(sTpsPbField);
  auto sTpsDivuView = make_view<double, 3>(sTpsDivuField);
  auto sQPbView = make_view<double, 3>(sQPbField);
  auto sQDivuView = make_view<double, 3>(sQDivuField);
  auto sQTpsuView = make_view<double, 3>(sQTpsuField);

  // Deserialize data
  for (size_t jw = 0; jw < trans.nw(); ++jw) {
    for (size_t jz2 = 0; jz2 < nflev; ++jz2) {
      for (size_t jz1 = 0; jz1 < nflev; ++jz1) {
        const size_t jv = jw*nflev*nflev + jz1*nflev + jz2;
        sDivPbView(jw, jz2, jz1) = sDivPb[jv];
        sQPbView(jw, jz2, jz1) = sQPb[jv];
        sQDivuView(jw, jz2, jz1) = sQDivu[jv];
      }
    }
    for (size_t jz2 = 0; jz2 < nflev+1; ++jz2) {
      for (size_t jz1 = 0; jz1 < nflev; ++jz1) {
        const size_t jv = jw*nflev*(nflev+1) + jz1*(nflev+1) + jz2;
        sTpsPbView(jw, jz2, jz1) = sTpsPb[jv];
        sTpsDivuView(jw, jz2, jz1) = sTpsDivu[jv];
      }
    }
    for (size_t jz2 = 0; jz2 < nflev; ++jz2) {
      for (size_t jz1 = 0; jz1 < nflev+1; ++jz1) {
        const size_t jv = jw*(nflev+1)*nflev + jz1*nflev + jz2;
        sQTpsuView(jw, jz2, jz1) = sQTpsu[jv];
      }
    }
  }

  oops::Log::trace() << "saber::bifourier::BifourierAromeLegacy::readBalance done" << std::endl;
}

// -----------------------------------------------------------------------------

void readCovariance(const eckit::mpi::Comm & comm,
                    const oops::Variables & activeVars,
                    const BifourierCovarianceReadParameters & params,
                    const BifourierTransform & trans,
                    atlas::FieldSet & fset) {
  oops::Log::trace() << "saber::bifourier::BifourierAromeLegacy::readCovariance starting"
    << std::endl;

  // Get number of levels
  const size_t nflev = activeVars["air_upward_absolute_vorticity"].getLevels();

  // Allocate local vectors
  std::vector<double> vorCov(trans.nw()*nflev*nflev);
  std::vector<double> divuCov(trans.nw()*nflev*nflev);
  std::vector<double> tPsuCov(trans.nw()*(nflev+1)*(nflev+1));
  std::vector<double> quCov(trans.nw()*nflev*nflev);

  // Tag root
  const size_t tagRoot = 123;

  if (comm.rank() == 0) {
    // Allocate global vectors
    std::vector<double> vorCovGlb(trans.nwGlb()*nflev*nflev);
    std::vector<double> divuCovGlb(trans.nwGlb()*nflev*nflev);
    std::vector<double> tPsuCovGlb(trans.nwGlb()*(nflev+1)*(nflev+1));
    std::vector<double> quCovGlb(trans.nwGlb()*nflev*nflev);

    if (params.inputFileFormat.value() == "arome legacy binary") {
      // Read Fortran unformatted file (from readjbdat96.F90)
      bifourier_arome_legacy_covariance_f90(params.toConfiguration(), trans.nwGlb(), nflev,
        vorCovGlb.data(), divuCovGlb.data(), tPsuCovGlb.data(), quCovGlb.data());
    } else if (params.inputFileFormat.value() == "arome legacy netcdf") {
      // NetCDF file path
      std::string ncFilePath = params.inputFile.value();

      // NetCDF IDs
      int ncid, retval, dimid, varid;
      size_t nflevFromFile, nwFromFile;

      // Open NetCDF file
      if ((retval = nc_open(ncFilePath.c_str(), NC_NOWRITE, &ncid))) ERR(retval, ncFilePath);

      // Check dimensions
      if ((retval = nc_inq_dimid(ncid, "NFLEV", &dimid))) ERR(retval, "NFLEV");
      if ((retval = nc_inq_dimlen(ncid, dimid, &nflevFromFile))) ERR(retval, "NFLEV");
      ASSERT(nflevFromFile == nflev);
      if ((retval = nc_inq_dimid(ncid, "NSMAXP1", &dimid))) ERR(retval, "NSMAXP1");
      if ((retval = nc_inq_dimlen(ncid, dimid, &nwFromFile))) ERR(retval, "NSMAXP1");
      ASSERT(nwFromFile == trans.nwGlb());

      // Get variables
      if ((retval = nc_inq_varid(ncid, "VOR_VERTCOV", &varid))) ERR(retval, "VOR_VERTCOV");
      if ((retval = nc_get_var_double(ncid, varid, vorCovGlb.data()))) ERR(retval, "VOR_VERTCOV");
      if ((retval = nc_inq_varid(ncid, "DIVU_VERTCOV", &varid))) ERR(retval, "DIVU_VERTCOV");
      if ((retval = nc_get_var_double(ncid, varid, divuCovGlb.data())))
        ERR(retval, "DIVU_VERTCOV");
      if ((retval = nc_inq_varid(ncid, "TPSU_VERTCOV", &varid))) ERR(retval, "TPSU_VERTCOV");
      if ((retval = nc_get_var_double(ncid, varid, tPsuCovGlb.data())))
        ERR(retval, "TPSU_VERTCOV");
      if ((retval = nc_inq_varid(ncid, "QU_VERTCOV", &varid))) ERR(retval, "QU_VERTCOV");
      if ((retval = nc_get_var_double(ncid, varid, quCovGlb.data()))) ERR(retval, "QU_VERTCOV");

      // Close file
      if ((retval = nc_close(ncid))) ERR(retval, ncFilePath);
    }

    // Print norms
    oops::Log::info() << "Info     : Covariance norms: " << std::endl;
    const double vorCovNorm = std::inner_product(vorCovGlb.begin(), vorCovGlb.end(),
      vorCovGlb.begin(), 0.0);
    oops::Log::test() << "- vorVertCov: " << vorCovNorm << std::endl;
    const double divuCovNorm = std::inner_product(divuCovGlb.begin(), divuCovGlb.end(),
      divuCovGlb.begin(), 0.0);
    oops::Log::test() << "- divuVertCov: " << divuCovNorm << std::endl;
    const double tPsuCovNorm = std::inner_product(tPsuCovGlb.begin(), tPsuCovGlb.end(),
      tPsuCovGlb.begin(), 0.0);
    oops::Log::test() << "- tPsuVertCov: " << tPsuCovNorm << std::endl;
    const double quCovNorm = std::inner_product(quCovGlb.begin(), quCovGlb.end(),
      quCovGlb.begin(), 0.0);
    oops::Log::test() << "- quVertCov: " << quCovNorm << std::endl;

    // Copy data
    for (size_t jw = 0; jw < trans.nw(); ++jw) {
      const size_t jwGlb = jw + trans.nwStart();
      for (size_t jz2 = 0; jz2 < nflev; ++jz2) {
        for (size_t jz1 = 0; jz1 < nflev; ++jz1) {
          const size_t index = jw*nflev*nflev + jz2*nflev + jz1;
          const size_t indexGlb = jwGlb*nflev*nflev + jz2*nflev + jz1;
          vorCov[index] = vorCovGlb[indexGlb];
          divuCov[index] = divuCovGlb[indexGlb];
          quCov[index] = quCovGlb[indexGlb];
        }
      }
      for (size_t jz2 = 0; jz2 < nflev+1; ++jz2) {
        for (size_t jz1 = 0; jz1 < nflev+1; ++jz1) {
          const size_t index = jw*(nflev+1)*(nflev+1) + jz2*(nflev+1) + jz1;
          const size_t indexGlb = jwGlb*(nflev+1)*(nflev+1) + jz2*(nflev+1) + jz1;
          tPsuCov[index] = tPsuCovGlb[indexGlb];
        }
      }
    }

    // Send data
    int tag = tagRoot;
    for (size_t jt = 1; jt < comm.size(); ++jt) {
      // Create vectors
      std::vector<double> vorCovSend(trans.nwPerTask()[jt]*nflev*nflev);
      std::vector<double> divuCovSend(trans.nwPerTask()[jt]*nflev*nflev);
      std::vector<double> tPsuCovSend(trans.nwPerTask()[jt]*(nflev+1)*(nflev+1));
      std::vector<double> quCovSend(trans.nwPerTask()[jt]*nflev*nflev);

      // Fill vectors
      for (size_t jwSend = 0; jwSend < trans.nwPerTask()[jt]; ++jwSend) {
        const size_t jwGlb = jwSend + trans.nwStartPerTask()[jt];
        for (size_t jz2 = 0; jz2 < nflev; ++jz2) {
          for (size_t jz1 = 0; jz1 < nflev; ++jz1) {
            const size_t indexSend = jwSend*nflev*nflev + jz2*nflev + jz1;
            const size_t indexGlb = jwGlb*nflev*nflev + jz2*nflev + jz1;
            vorCovSend[indexSend] = vorCovGlb[indexGlb];
            divuCovSend[indexSend] = divuCovGlb[indexGlb];
            quCovSend[indexSend] = quCovGlb[indexGlb];
          }
        }
      }
      for (size_t jwSend = 0; jwSend < trans.nwPerTask()[jt]; ++jwSend) {
        const size_t jwGlb = jwSend + trans.nwStartPerTask()[jt];
        for (size_t jz2 = 0; jz2 < nflev+1; ++jz2) {
          for (size_t jz1 = 0; jz1 < nflev+1; ++jz1) {
            const size_t indexSend = jwSend*(nflev+1)*(nflev+1) + jz2*(nflev+1) + jz1;
            const size_t indexGlb = jwGlb*(nflev+1)*(nflev+1) + jz2*(nflev+1) + jz1;
            tPsuCovSend[indexSend] = tPsuCovGlb[indexGlb];
          }
        }
      }

      // Send vectors
      comm.send(vorCovSend.data(), trans.nwPerTask()[jt]*nflev*nflev, jt, tag+0);
      comm.send(divuCovSend.data(), trans.nwPerTask()[jt]*nflev*nflev, jt, tag+1);
      comm.send(tPsuCovSend.data(), trans.nwPerTask()[jt]*(nflev+1)*(nflev+1), jt, tag+2);
      comm.send(quCovSend.data(), trans.nwPerTask()[jt]*nflev*nflev, jt, tag+3);

      // Update tag
      tag += 4;
    }
  }

  // Receive vectors
  if (comm.rank() > 0) {
    const int tagBase = tagRoot+4*(comm.rank()-1);
    comm.receive(vorCov.data(), trans.nw()*nflev*nflev, 0, tagBase+0);
    comm.receive(divuCov.data(), trans.nw()*nflev*nflev, 0, tagBase+1);
    comm.receive(tPsuCov.data(), trans.nw()*(nflev+1)*(nflev+1), 0, tagBase+2);
    comm.receive(quCov.data(), trans.nw()*nflev*nflev, 0, tagBase+3);
  }

  // MPI barrier
  comm.barrier();

  // Constant coefficient
  const double zmovern = static_cast<double>(trans.ellips().size())
    / static_cast<double>(trans.nwGlb()-1);

  // Theoretical REDNMC value
  const double sqrtHalf = std::sqrt(0.5);

  for (const auto & var : activeVars) {
    // Get number of levels
    const size_t nz = var.getLevels();

    // Create covariance field
    createField3D("cov", trans.nw(), var, fset);

    // Get covariance view
    auto covView = getView3D("cov", var, fset);

    // Copy data
    for (size_t jw = 0; jw < trans.nw(); ++jw) {
      for (size_t jz2 = 0; jz2 < nz; ++jz2) {
        for (size_t jz1 = 0; jz1 < nz; ++jz1) {
          const size_t index = jw*nz*nz + jz2*nz + jz1;
          if (var.name() == "air_upward_absolute_vorticity") {
            covView(jw, jz2, jz1) = vorCov[index];
          }
          if (var.name() == "air_horizontal_divergence") {
            covView(jw, jz2, jz1) = divuCov[index];
          }
          if (var.name() == "air_temperature_and_log_of_air_pressure_at_surface") {
            covView(jw, jz2, jz1) = tPsuCov[index];
          }
          if (var.name() == "water_vapor_mixing_ratio_wrt_moist_air") {
            covView(jw, jz2, jz1) = quCov[index];
          }
        }
      }
    }

    // Get correlation square-root view
    auto corSqrtView = getView3D("corSqrt", var, fset);

    // From suejbstd.F90

    // Create standard-deviation profiles
    std::vector<double> vertStd(nz, 0.0);

    // Compute vertical variance for each level
    for (size_t jw = 0; jw < trans.nwRoot(); ++jw) {
      for (size_t jz = 0; jz < nz; ++jz) {
        vertStd[jz] += covView(jw, jz, jz);
      }
    }

    // Communication
    comm.allReduceInPlace(vertStd.begin(), vertStd.end(), eckit::mpi::sum());

    // Take variance square-root
    for (size_t jz = 0; jz < nz; ++jz) {
      vertStd[jz] = std::sqrt(vertStd[jz]);
    }

    // From suejbcor.F90

    // Compute vertical correlation square-root
    const double zeps = 1.0e-99;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    for (size_t jw = 0; jw < trans.nw(); ++jw) {
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

    // Set wavenumber 0 to zero for vorticity and divergence
    if (var.name() == "air_upward_absolute_vorticity"
      || var.name() == "air_horizontal_divergence") {
      for (size_t jw = 0; jw < trans.nw(); ++jw) {
        const size_t jwGlb = jw + trans.nwStart();
        if (jwGlb == 0) {
          for (size_t jz2 = 0; jz2 < nz; ++jz2) {
            for (size_t jz1 = 0; jz1 < nz; ++jz1) {
              covView(jw, jz2, jz1) = 0.0;
            }
          }
        }
      }
    }

    // Compute sum over wavenumbers
    std::vector<double> sum(nz, 0.0);
    for (size_t jw = 0; jw < trans.nwRoot(); ++jw) {
      for (size_t jz = 0; jz < nz; ++jz) {
        sum[jz] += covView(jw, jz, jz);
      }
    }

    // Communication
    comm.allReduceInPlace(sum.begin(), sum.end(), eckit::mpi::sum());

    // Create horizontal covariance field
    atlas::Field horCovField("horCov", make_datatype<double>(),
      make_shape(trans.nw(), nz));

    // Get horizontal covariance view
    auto horCovView = make_view<double, 2>(horCovField);

    // Normalize sums
    for (size_t jw = 0; jw < trans.nw(); ++jw) {
      for (size_t jz = 0; jz < nz; ++jz) {
        horCovView(jw, jz) = covView(jw, jz, jz)/sum[jz];
      }
    }

    // Apply weight
    double zWeight;
    for (size_t jw = 0; jw < trans.nw(); ++jw) {
      const size_t jwGlb = jw + trans.nwStart();
      if (jwGlb != 0 && jwGlb != trans.nwGlb()-1) {
        zWeight = 2.0*M_PI*static_cast<double>(jwGlb)*zmovern;
      } else if (jwGlb == 0) {
        zWeight = M_PI*zmovern/4.0;
      } else if (jwGlb == trans.nwGlb()-1) {
        zWeight = M_PI*(static_cast<double>(trans.nwGlb()-1)-0.25)*zmovern;
      }
      for (size_t jz = 0; jz < nz; ++jz) {
        horCovView(jw, jz) = std::max(1.0e-20, std::sqrt(horCovView(jw, jz) / zWeight));
      }
    }

    // Merge contributions
    for (size_t jw = 0; jw < trans.nw(); ++jw) {
      for (size_t jz2 = 0; jz2 < nz; ++jz2) {
        for (size_t jz1 = 0; jz1 < nz; ++jz1) {
          corSqrtView(jw, jz2, jz1) *= horCovView(jw, jz2) * vertStd[jz2] * sqrtHalf;
        }
      }
    }

    // Get standard-deviation view
    auto stdDevView = getViewProfile("stdDev", var, fset);

    // Compute variance
    std::vector<double> variance(nz, 0.0);
    for (size_t jz2 = 0; jz2 < nz; ++jz2) {
      for (size_t js = 0; js < trans.ns(); ++js) {
        if (trans.jq(js) == 0) {
          const size_t jw = trans.jw(js);
          for (size_t jz1 = 0; jz1 < nz; ++jz1) {
            variance[jz2] += corSqrtView(jw, jz2, jz1)*corSqrtView(jw, jz2, jz1)
              *trans.spNorm(js)*trans.spNorm(js);
          }
        }
      }
    }

    // Communication
    comm.allReduceInPlace(variance.begin(), variance.end(), eckit::mpi::sum());

    // Compute standard-deviation
    for (size_t jz = 0; jz < nz; ++jz) {
      stdDevView(jz) = std::sqrt(variance[jz]);
    }

    // Apply inverse standard-deviation to normalize correlation square-root
    for (size_t jz2 = 0; jz2 < nz; ++jz2) {
      ASSERT(stdDevView(jz2) > 0.0);
      const double norm = 1.0/stdDevView(jz2);
      for (size_t jw = 0; jw < trans.nw(); ++jw) {
        for (size_t jz1 = 0; jz1 < nz; ++jz1) {
          corSqrtView(jw, jz2, jz1) *= norm;
        }
      }
    }
  }

  oops::Log::trace() << "saber::bifourier::BifourierAromeLegacy::readCovariance done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace arome_legacy
}  // namespace bifourier
}  // namespace saber

