/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <memory>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/grid.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Timer.h"

#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/oops/SaberOuterBlockParametersBase.h"
#include "saber/spectralb/gaussutils.h"
#include "saber/spectralb/SPTOUV.h"

namespace saber {
namespace spectralb {

namespace {

SPTOUVParameters createSPTOUVParams(const eckit::Configuration & conf) {
  SPTOUVParameters params;
  params.deserialize(conf);
  return params;
}

oops::Variables createInputVars(const SPTOUVParameters & params) {
  oops::Variables inputVars;
  oops::Variables outputVars = params.outputVars.value();

  for (auto & var : outputVars.variables()) {
    if (var.compare("eastward_wind") == 0 || var.compare("northward_wind") == 0) {
    } else {
      inputVars.push_back(var);
    }
  }
  inputVars.push_back("streamfunction");
  inputVars.push_back("velocity_potential");
  //inputVars.push_back("vorticity");
  //inputVars.push_back("divergence");

  return inputVars;
}


oops::Variables createActiveVars(const oops::Variables & inputVars,
                                 const oops::Variables & outputVars) {
  oops::Variables activeVars;
  for (auto & var : outputVars.variables()) {
    if (!inputVars.has(var)){
      activeVars.push_back(var);
    }
  }
  for (auto & var : inputVars.variables()) {
    if (!outputVars.has(var) && !activeVars.has(var)){
      activeVars.push_back(var);
    }
  }

  return activeVars;
}


double norm(const atlas::FieldSet & fset) {
  double x(0.0);
  for (auto & f : fset) {
     auto fldView = atlas::array::make_view<double, 2>(f);
     for (atlas::idx_t i = 0; i < f.shape(0)   ; ++i) {
       for (atlas::idx_t j = 0; j < f.shape(1)   ; ++j) {
         x += fldView(i, j) *  fldView(i, j);
       }

     }
  }
  return x;
}

}
// -----------------------------------------------------------------------------

static SaberOuterBlockMaker<SPTOUV> makerSPTOUV_("SPTOUV");

// Build input functionspace from output functionspace
// It is the output functionspace that is in the argument.
// So we need to create spectral functionspace !!!

// We need to write separate functions that allow

// Rule:
// Expect each saber block in its multiply method as a final step to apply a
// halo exchange at the end of the method if it needs it.
// and at the beginning of the adjoint of the saber block to apply the adjoint
// halo exchange.

// "active variables" - now required in yaml
// "input variables" - optional in yaml - input for the multiply
//                   - sum active and passive variables
SPTOUV::SPTOUV(const eckit::mpi::Comm & comm,
               const atlas::FunctionSpace & outputFunctionSpace,
               const atlas::FieldSet & outputExtraFields,
               const std::vector<size_t> & activeVariableSizes,
               const eckit::Configuration & conf,
               const atlas::FieldSet & xb,
               const atlas::FieldSet & fg,
               const std::vector<atlas::FieldSet> & fsetVec)
  : SaberOuterBlockBase(conf),
    params_(createSPTOUVParams(conf)),
    outputFunctionSpace_(atlas::FunctionSpace(outputFunctionSpace)),
    outputVars_(params_.outputVars.value()),
    inputVars_(createInputVars(params_)),
    activeVariableSizes_(activeVariableSizes),
    gaussGrid_(atlas::StructuredGrid(params_.gaussGridUid)),
    specFS_(createSpectralFunctionSpace(gaussGrid_, activeVariableSizes)),
    inputFunctionSpace_(atlas::FunctionSpace(specFS_)),
    transFS_(atlas::trans::Trans(outputFunctionSpace_, inputFunctionSpace_))

{
  oops::Log::trace() << classname() << "::SPTOUV starting" << std::endl;
  SaberOuterBlockBase::inputFunctionSpace_ = inputFunctionSpace_;
  SaberOuterBlockBase::inputVars_ = inputVars_;
  SaberOuterBlockBase::inputExtraFields_ = outputExtraFields;

  std::cout << "SPTOUV inputVars " << inputVars_.variables() << std::endl;

  oops::Log::trace() << classname() << "::SPTOUV done" << std::endl;
}

// -----------------------------------------------------------------------------
SPTOUV::~SPTOUV() {
  oops::Log::trace() << classname() << "::~SPTOUV starting" << std::endl;
  util::Timer timer(classname(), "~SPTOUV");
  oops::Log::trace() << classname() << "::~SPTOUV done" << std::endl;
}

// -----------------------------------------------------------------------------
void SPTOUV::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;

  std::cout << "inputVars multiply start " <<
               inputVars_.variables() << " " <<  norm(fset) << std::endl;

  for (atlas::idx_t jvar = 0; jvar < fset.size(); ++jvar) {
    auto view = atlas::array::make_view<double, 2>(fset[jvar]);
    for (atlas::idx_t jnode = 0; jnode < fset[jvar].shape(0); ++jnode) {
      std::cout << "spec start mult = " << fset[jvar].name()
                << " " << jnode << " " << view(jnode, 0) << std::endl;
    }
  }

  // Create empty Model fieldset
  atlas::FieldSet newFields = atlas::FieldSet();

  // copy "passive variables"
  std::vector<std::string> fsetNames = fset.field_names();

  std::cout << "fset fieldnames" << fsetNames << std::endl;

  for (auto & s : fset.field_names()) {
     if (outputVars_.has(s)) {
       newFields.add(fset[s]);
     }
  }

  // check variables and if streamfunction, velocity potential scale
  // by n(n+1) / a   and rename fields to vorticity_spectral_2D
  // divergence_spectral_2D
  const int N = atlas::GaussianGrid(gaussGrid_).N();

  if (inputVars_.has("streamfunction") && inputVars_.has("velocity_potential")) {
    applyNtimesNplus1SpectralScaling(
      oops::Variables(std::vector<std::string>({"streamfunction", "velocity_potential"})),
      oops::Variables(std::vector<std::string>({"vorticity", "divergence"})),
      specFS_, N, fset);
  }

  std::cout << activeVariableSizes_.size() << std::endl;
  std::cout << outputFunctionSpace_.type() << std::endl;
  std::cout << fset.field_names() << std::endl;

  atlas::Field uvgp = allocateGaussUVField(outputFunctionSpace_,
                                           inputVars_, activeVariableSizes_);

  std::cout << "after allocate UV  " << fset.field_names() << std::endl;
  for (atlas::Field & f  : fset) {
    std::cout << f.name() << std::endl;
  }

  std::cout << fset["divergence"].name() << std::endl;

  // transform to gaussian grid
  transFS_.invtrans_vordiv2wind(fset["vorticity"], fset["divergence"], uvgp);

  atlas::FieldSet uvfset = convertUVToFieldSet(uvgp);

  newFields.add(uvfset["eastward_wind"]);
  newFields.add(uvfset["northward_wind"]);

  std::cout << "after nf" << std::endl;

  fset = newFields;

  for (atlas::idx_t jvar = 0; jvar < fset.size(); ++jvar) {
    auto view = atlas::array::make_view<double, 2>(fset[jvar]);
    for (atlas::idx_t jnode = 0; jnode < fset[jvar].shape(0); ++jnode) {
      std::cout << "multiply end  = " << fset[jvar].name()
                << " " << jnode << " " << view(jnode, 0) << std::endl;
    }
  }


  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}


// -----------------------------------------------------------------------------
void SPTOUV::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << " "  << norm(fset) << std::endl;

  for (atlas::idx_t jvar = 0; jvar < fset.size(); ++jvar) {
    auto view = atlas::array::make_view<double, 2>(fset[jvar]);
    for (atlas::idx_t jnode = 0; jnode < fset[jvar].shape(0); ++jnode) {
      std::cout << "multiplyAD start  = " << fset[jvar].name()
                << " " << jnode << " " << view(jnode, 0) << std::endl;
    }
  }

  // Create empty Model fieldset
  atlas::FieldSet newFields = atlas::FieldSet();

  // copy "passive variables"
  std::vector<std::string> fsetNames = fset.field_names();

  for (auto & s : fset.field_names()) {
     if (inputVars_.has(s)) {
       newFields.add(fset[s]);
     }
  }

  std::cout << "adj fn space type newFields.size" << fset[0].functionspace().type() << " " <<
               newFields.size()  << std::endl;

  atlas::Field uvgp = convertUVToFieldSetAD(fset);

  atlas::FieldSet specfset = allocateSpectralVortDiv(specFS_, inputVars_, activeVariableSizes_);

  transFS_.invtrans_vordiv2wind_adj(uvgp, specfset["vorticity"], specfset["divergence"]);

  const int N = atlas::GaussianGrid(gaussGrid_).N();

  if (inputVars_.has("streamfunction") && inputVars_.has("velocity_potential")) {
    applyNtimesNplus1SpectralScaling(
      oops::Variables(std::vector<std::string>({"vorticity", "divergence"})),
      oops::Variables(std::vector<std::string>({"streamfunction", "velocity_potential"})),
      specFS_, N, specfset);
  }

  newFields.add(specfset[0]);
  newFields.add(specfset[1]);

  fset = newFields;

  for (atlas::idx_t jvar = 0; jvar < fset.size(); ++jvar) {
    auto view = atlas::array::make_view<double, 2>(fset[jvar]);
    for (atlas::idx_t jnode = 0; jnode < fset[jvar].shape(0); ++jnode) {
      std::cout << "multiplyAD done = " << fset[jvar].name()
                << " " << jnode << " " << view(jnode, 0) << std::endl;
    }
  }

  oops::Log::trace() << classname() << "::multiplyAD done" << " "  << norm(fset) << std::endl;
}

// -----------------------------------------------------------------------------
void SPTOUV::calibrationInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::calibrationInverseMultiply starting" << std::endl;
  oops::Log::trace() << classname() << "::calibrationInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------
void SPTOUV::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace spectralb
}  // namespace saber
