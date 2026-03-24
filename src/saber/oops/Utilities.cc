/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "oops/base/FieldSet3D.h"
#include "oops/util/for_each.h"

#include "saber/oops/Utilities.h"

namespace saber {

// -----------------------------------------------------------------------------

// Return inner variables as outer variables + inner active variables
// Can be used to help define innerVars_ member in SABER outer blocks
oops::Variables getUnionOfInnerActiveAndOuterVars(const SaberBlockParametersBase & params,
                                                  const oops::Variables & outerVars) {
  oops::Variables innerVars(outerVars);
  innerVars += params.activeInnerVars(outerVars);
  return innerVars;
}

// -----------------------------------------------------------------------------

// Return inner variables that are not outer variables
oops::Variables getInnerOnlyVars(const SaberBlockParametersBase & params,
                                 const oops::Variables & outerVars) {
  oops::Variables innerOnlyVars(getUnionOfInnerActiveAndOuterVars(params, outerVars));
  innerOnlyVars -= outerVars;
  return innerOnlyVars;
}

// -----------------------------------------------------------------------------

void setMPI(eckit::LocalConfiguration & conf,
            const int & mpi) {
  oops::Log::trace() << "setMPI starting" << std::endl;

  if (conf.has("mpi pattern")) {
    std::string mpiPattern = conf.getString("mpi pattern");
    util::seekAndReplace(conf, mpiPattern, std::to_string(mpi));
  }

  oops::Log::trace() << "setMPI done" << std::endl;
}

// -----------------------------------------------------------------------------

void checkFieldsAreNotAllocated(const oops::FieldSet3D & fset,
                                const oops::Variables & vars) {
  for (const auto& var : vars) {
    if (fset.has(var.name())) {
      throw eckit::UserError("Variable " + var.name() + " is already allocated in FieldSet.",
                             Here());
    }
  }
}

// -----------------------------------------------------------------------------

void allocateMissingFields(oops::FieldSet3D & fset,
                           const oops::Variables & varsToAllocate,
                           const oops::Variables & varsWithLevels,
                           const atlas::FunctionSpace & functionSpace) {
  oops::Log::trace() << "allocateMissingFields starting" << std::endl;
  for (const auto& var : varsToAllocate) {
    if (!fset.has(var.name())) {
      oops::Log::info() << "Info     : Allocating " << var.name() << std::endl;
      auto field = functionSpace.createField<double>(
                atlas::option::name(var.name()) |
                atlas::option::levels(varsWithLevels[var.name()].getLevels()));
      atlas::array::make_view<double, 2>(field).assign(0.0);
      field.set_dirty(false);
      fset.add(field);
    }
  }
  oops::Log::trace() << "allocateMissingFields done" << std::endl;
}

// -----------------------------------------------------------------------------

size_t getNensFromConfig(const eckit::Configuration & conf) {
  size_t nens = 0;
  for (const auto & ensType : {"ensemble", "ensemble pert", "ensemble base",
    "ensemble pert on other geometry"}) {
    if (conf.has(ensType)) {
      eckit::LocalConfiguration ensTypeConf = conf.getSubConfiguration(ensType);

      ASSERT(ensTypeConf.has("members from template") || ensTypeConf.has("members"));
      ASSERT(!(ensTypeConf.has("members from template") && ensTypeConf.has("members")));
      nens = getNensFromConfig(ensTypeConf);
    }
  }
  if (conf.has("members")) {
    const auto members = conf.getSubConfigurations("members");
    nens = members.size();
  } else if (conf.has("members from template")) {
    const auto members = conf.getSubConfiguration("members from template");
    ASSERT(members.has("nmembers"));
    ASSERT(members.has("pattern"));
    ASSERT(members.has("template"));
    nens = members.getInt("nmembers");
  }
  return nens;
}

// -----------------------------------------------------------------------------

eckit::LocalConfiguration getEnsSubconfig(const eckit::Configuration & conf, size_t iens) {
  // expecting either `members` (list) or `members from template` (object with nmembers/template).
  eckit::LocalConfiguration memConf;
  if (conf.has("members")) {
    const auto members = conf.getSubConfigurations("members");
    if (!members.empty()) memConf = members[iens];
  } else {
    eckit::LocalConfiguration members = conf.getSubConfiguration("members from template");
    memConf = members.getSubConfiguration("template");
    const std::string pattern = members.getString("pattern");
    const int zpad = members.getInt("zero padding", 0);
    const std::vector<size_t> except = members.getUnsignedVector("except", {});
    size_t index = members.getUnsigned("start", 1);
    for (size_t jj = 0; jj <= iens; ++jj) {
      // Check for excluded members
      while (std::count(except.begin(), except.end(), index)) {
        index++;
      }
      // Update counter
      if (jj < iens) index++;
    }
    util::seekAndReplace(memConf, pattern, index, zpad);
  }
  return memConf;
}

// -----------------------------------------------------------------------------

void cvToFset(const atlas::Field & cv,
              oops::FieldSet3D & fset,
              const size_t & offset,
              const oops::Variables & vars) {
  // Copy from control vector to fieldset
  size_t fld_offset = offset;
  const auto cvView = atlas::array::make_view<double, 1>(cv);
  for (const auto & var : vars) {
    auto & field = fset[var.name()];
    size_t vt_nodes = field.shape(1);
    auto view = atlas::array::make_view<double, 2>(field);
    util::for_each_index(
      util::make_index_space_2d(util::IndexRange::include_halo, field),
      [=](atlas::idx_t hz_index, atlas::idx_t vt_index) mutable {
        size_t index = fld_offset + vt_index + hz_index*vt_nodes;
        view(hz_index, vt_index) = cvView(index);
      });
    fld_offset += field.size();
  }
}

// -----------------------------------------------------------------------------

void fsetToCv(const oops::FieldSet3D & fset,
              atlas::Field & cv,
              const size_t & offset,
              const oops::Variables & vars) {
  // Copy from fieldset to control vector
  size_t fld_offset = offset;
  auto cvView = atlas::array::make_view<double, 1>(cv);
  for (const auto & var : vars) {
    const auto & field = fset[var.name()];
    const auto view = atlas::array::make_view<double, 2>(field);
    size_t vt_nodes = field.shape(1);
    util::for_each_index(
      util::make_index_space_2d(util::IndexRange::include_halo, field),
      [=](atlas::idx_t hz_index, atlas::idx_t vt_index) mutable {
        size_t index = fld_offset + vt_index + hz_index*vt_nodes;
        cvView(index) = view(hz_index, vt_index);
      });
    fld_offset += field.size();
  }
}

// -----------------------------------------------------------------------------

}  // namespace saber
