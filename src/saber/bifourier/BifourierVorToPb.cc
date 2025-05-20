/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#include "saber/bifourier/BifourierVorToPb.h"

#include <netcdf.h>

#include <algorithm>

#include "atlas/field.h"

#include "saber/bifourier/BifourierAromeLegacy.h"

#define ERR(e, msg) {std::string s(nc_strerror(e)); \
  throw eckit::Exception(s + " : " + msg, Here());}

using atlas::array::make_view;

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<BifourierVorToPb> makerBifourierVorToPb_("BifourierVorToPb");

// -----------------------------------------------------------------------------

static std::vector<double> fact1Static;

// -----------------------------------------------------------------------------

std::vector<double> & BifourierVorToPb::fact1() {
  return fact1Static;
}

// -----------------------------------------------------------------------------

BifourierVorToPb::BifourierVorToPb(const oops::GeometryData & outerGeometryData,
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
  oops::Log::trace() << classname() << "::BifourierVorToPb starting" << std::endl;

  // Inner variables
  if (!params_.backward.value()) {
    // Forward mode
    innerVars_ -= innerVars_["balanced_air_pressure"];
  } else {
    // Backward mode
    innerVars_.push_back("balanced_air_pressure");
    innerVars_["balanced_air_pressure"].setLevels(
      innerVars_["air_upward_absolute_vorticity"].getLevels());
  }

  // Retrieve spectral transform
  trans_ = transStore_.retrieveTransform(outerGeometryData);

  if (params_.read.value() == boost::none && fact1().size() == 0) {
    // Allocate fact1
    fact1().resize(trans_->ns());

    // Get change of variable parameters from configuration or from spectral transform
    const auto & nkFromConf = params_.nk.value();
    const size_t nk = (nkFromConf != boost::none) ? *nkFromConf : trans_->nk();
    const auto & nlFromConf = params_.nl.value();
    const size_t nl = (nlFromConf != boost::none) ? *nlFromConf : trans_->nl();
    const auto & meanLatFromConf = params_.meanLat.value();
    const double meanLat = (meanLatFromConf != boost::none) ? *meanLatFromConf : trans_->meanLat();

    // Compute change of variable factor
    const size_t nwGlb = std::max(nk-1, nl-1);
    const double zromega = 0.7292115e-4;
    const double zcc = -2.0*zromega*std::sin(meanLat*M_PI/180.0);
    const double zly = 2.0*static_cast<double>(nwGlb)*trans_->dy();
    const double zfact1 = zcc*(zly/(2.0*M_PI))*(zly/(2.0*M_PI));
    for (size_t js = 0; js < trans_->ns(); ++js) {
      const double kstar = trans_->kstar(trans_->jk(js), trans_->jl(js), nk, nl, nwGlb);
      if (kstar > 0.0) {
        fact1()[js] = zfact1/(kstar*kstar);
      } else {
        fact1()[js] = 0.0;
      }
    }
  }

  oops::Log::trace() << classname() << "::BifourierVorToPb done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierVorToPb::multiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  if (!params_.backward.value()) {
    // Get inner field
    const auto vorField = fset["air_upward_absolute_vorticity"];
    const size_t nz = vorField.levels();

    // Create outer field
    atlas::Field pbField = trans_->spFspace()->createField<double>(
      atlas::option::name("balanced_air_pressure") | atlas::option::levels(nz));

    // Get fields views
    const auto vorView = make_view<double, 2>(vorField);
    auto pbView = make_view<double, 2>(pbField);

    // Apply change of variable
    for (size_t js = 0; js < trans_->ns(); ++js) {
      for (size_t jz = 0; jz < nz; ++jz) {
        pbView(js, jz) = vorView(js, jz)*fact1()[js];
      }
    }

    // Add outer field
    fset.add(pbField);
  } else {
    // Remove outer field
    util::removeFieldsFromFieldSet(fset.fieldSet(), {"balanced_air_pressure"});
  }

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierVorToPb::multiplyAD(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  if (!params_.backward.value()) {
    // Get fields
    const auto pbField = fset["balanced_air_pressure"];
    auto vorField = fset["air_upward_absolute_vorticity"];
    const size_t nz = vorField.levels();

    // Get fields views
    const auto pbView = make_view<double, 2>(pbField);
    auto vorView = make_view<double, 2>(vorField);

    // Apply change of variable, adjoint
    for (size_t js = 0; js < trans_->ns(); ++js) {
      for (size_t jz = 0; jz < nz; ++jz) {
        vorView(js, jz) += pbView(js, jz)*fact1()[js];
      }
    }

    // Remove outer field
    util::removeFieldsFromFieldSet(fset.fieldSet(), {"balanced_air_pressure"});
  } else {
    // Get outer field
    const auto vorField = fset["air_upward_absolute_vorticity"];
    const size_t nz = vorField.levels();

    // Create inner field
    atlas::Field pbField = trans_->spFspace()->createField<double>(
      atlas::option::name("balanced_air_pressure") | atlas::option::levels(nz));

    // Get inner field view
    auto pbView = make_view<double, 2>(pbField);

    // Set inner field to zero
    pbView.assign(0.0);

    // Add outer field
    fset.add(pbField);
  }

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierVorToPb::leftInverseMultiply(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::leftInverseMultiply starting" << std::endl;

  if (!params_.backward.value()) {
    // Get fields
    const auto pbField = fset["balanced_air_pressure"];
    auto vorField = fset["air_upward_absolute_vorticity"];
    const size_t nz = vorField.levels();

    // Get fields views
    const auto pbView = make_view<double, 2>(pbField);
    auto vorView = make_view<double, 2>(vorField);

    // Apply change of variable, inverse
    for (size_t js = 0; js < trans_->ns(); ++js) {
      for (size_t jz = 0; jz < nz; ++jz) {
        if (std::abs(fact1()[js]) > 0.0) {
          vorView(js, jz) = pbView(js, jz)/fact1()[js];
        }
      }
    }

    // Remove outer field
    util::removeFieldsFromFieldSet(fset.fieldSet(), {"balanced_air_pressure"});
  } else {
    // Get inner field
    const auto vorField = fset["air_upward_absolute_vorticity"];
    const size_t nz = vorField.levels();

    // Create outer field
    atlas::Field pbField = trans_->spFspace()->createField<double>(
      atlas::option::name("balanced_air_pressure") | atlas::option::levels(nz));

    // Get fields views
    const auto vorView = make_view<double, 2>(vorField);
    auto pbView = make_view<double, 2>(pbField);

    // Apply change of variable
    for (size_t js = 0; js < trans_->ns(); ++js) {
      for (size_t jz = 0; jz < nz; ++jz) {
        pbView(js, jz) = vorView(js, jz)*fact1()[js];
      }
    }

    // Add outer field
    fset.add(pbField);
  }

  oops::Log::trace() << classname() << "::leftInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierVorToPb::read() {
  oops::Log::trace() << classname() << "::read starting" << std::endl;

  if (params_.read.value()->inputFile.value() != boost::none) {
    // Allocate fact1
    fact1().resize(trans_->ns());

    // Read data
    if (params_.read.value()->inputFileFormat.value() == "arome legacy binary"
      || params_.read.value()->inputFileFormat.value() == "arome legacy netcdf") {
      // Read fact1 from file
      arome_legacy::readVorToPb(comm_, *params_.read.value(), *trans_, fact1());
    } else {
      // NetCDF file path
      std::string ncFilePath = *params_.read.value()->inputFile.value();

      // NetCDF IDs
      int ncid, retval, nsGlb_id, varid;
      size_t nsGlbFromFile;

      // Allocate global vector
      std::vector<double> fact1Glb;
      if (comm_.rank() == 0) {
        fact1Glb.resize(trans_->nsGlb());
      }

      if (comm_.rank() == 0) {
        // Open NetCDF file
        if ((retval = nc_open(ncFilePath.c_str(), NC_NOWRITE, &ncid))) ERR(retval, ncFilePath);

        // Check dimension
        if ((retval = nc_inq_dimid(ncid, "nsGlb", &nsGlb_id))) ERR(retval, "nsGlb");
        if ((retval = nc_inq_dimlen(ncid, nsGlb_id, &nsGlbFromFile))) ERR(retval, "nsGlb");
        ASSERT(nsGlbFromFile == trans_->nsGlb());

        // Get variable ID
        if ((retval = nc_inq_varid(ncid, "fact1", &varid))) ERR(retval, "fact1");

        // Read data
        std::vector<double> fact1GlbOrdered(trans_->nsGlb());
        if ((retval = nc_get_var_double(ncid, varid, fact1GlbOrdered.data()))) ERR(retval, "fact1");

        // Reorder data
        for (size_t jsGlb = 0; jsGlb < trans_->nsGlb(); ++jsGlb) {
          fact1Glb[jsGlb] = fact1GlbOrdered[trans_->sMapping()[jsGlb]];
        }

        // Close file
        if ((retval = nc_close(ncid))) ERR(retval, ncFilePath);
      }

      // Scatter vector
      comm_.scatterv(fact1Glb.cbegin(), fact1Glb.cend(), trans_->sCounts(), trans_->sDispls(),
        fact1().begin(), fact1().end(), 0);
    }
  }

  oops::Log::trace() << classname() << "::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierVorToPb::write() const {
  oops::Log::trace() << classname() << "::write starting" << std::endl;

  if (params_.outputFile.value() != boost::none) {
    // NetCDF IDs
    int retval, ncid, nsGlb_id, d1D_id[1], var_id;

    // NetCDF file path
    std::string ncFilePath = *params_.outputFile.value();

    // Allocate global vector
    std::vector<double> fact1Glb;
    if (comm_.rank() == 0) {
      fact1Glb.resize(trans_->nsGlb());
    }

    // Gather data
    comm_.gatherv(fact1().cbegin(), fact1().cend(), fact1Glb.begin(), fact1Glb.end(),
      trans_->sCounts(), trans_->sDispls(), 0);

    if (comm_.rank() == 0) {
      // Create NetCDF file
      if ((retval = nc_create(ncFilePath.c_str(), NC_64BIT_OFFSET | NC_CLOBBER, &ncid)))
        ERR(retval, ncFilePath);

      // Create dimension
      if ((retval = nc_def_dim(ncid, "nsGlb", trans_->nsGlb(), &nsGlb_id))) ERR(retval, "nsGlb");

      // Dimensions array
      d1D_id[0] = nsGlb_id;

      // Define variable
      if ((retval = nc_def_var(ncid, "fact1", NC_DOUBLE, 1, d1D_id, &var_id))) ERR(retval, "fact1");

      // End definition mode
      if ((retval = nc_enddef(ncid))) ERR(retval, ncFilePath);

      // Reorder data
      std::vector<double> fact1GlbOrdered(trans_->nsGlb());
      for (size_t jsGlb = 0; jsGlb < trans_->nsGlb(); ++jsGlb) {
        fact1GlbOrdered[trans_->sMapping()[jsGlb]] = fact1Glb[jsGlb];
      }

      // Write data
      if ((retval = nc_put_var_double(ncid, var_id, fact1GlbOrdered.data()))) ERR(retval, "fact1");

      // Close file
      if ((retval = nc_close(ncid))) ERR(retval, ncFilePath);
    }
  }

  oops::Log::trace() << classname() << "::write done" << std::endl;
}

// -----------------------------------------------------------------------------

void BifourierVorToPb::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
