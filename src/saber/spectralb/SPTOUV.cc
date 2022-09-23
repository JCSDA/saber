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
#include "oops/util/Timer.h"
#include "oops/util/abor1_cpp.h"

#include "saber/oops/SaberOuterBlockBase.h"
#include "saber/oops/SaberOuterBlockParametersBase.h"
#include "saber/spectralb/gaussutils.h"
#include "saber/spectralb/SPTOUV.h"

namespace saber {
namespace spectralb {
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
    outputFunctionSpace_(outputFunctionSpace),
    gaussGrid_(),
    transFS_()
{
  oops::Log::trace() << classname() << "::SPTOUV starting" << std::endl;

  // Deserialize configuration
  SPTOUVParameters params;
  params.deserialize(conf);

  inputVars_ = params.outputVars.value();
  specFS_ = createSpectralFunctionSpace(gaussGrid_, activeVariableSizes);
  inputFunctionSpace_ = atlas::FunctionSpace(specFS_);
  inputExtraFields_ = outputExtraFields;

  gaussGrid_ = atlas::StructuredGrid(params.gaussGridUid);
  activeVars_ = params.activeVars.value().get();
  transFS_ = atlas::trans::Trans(outputFunctionSpace_, inputFunctionSpace_);


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

  // Create empty Model fieldset
  atlas::FieldSet newFields = atlas::FieldSet();

  // copy "passive variables"
  std::vector<std::string> fsetNames = fset.field_names();
  for (auto s : fset.field_names()) {
    if (!activeVars_.has(s)) {
      newFields.add(fset[s]);
    }
  }

  // check variables and if streamfunction, velocity potential scale
  // by n(n+1) / a   and rename fields to vorticity_spectral_2D
  // divergence_spectral_2D
  const int N = atlas::GaussianGrid(gaussGrid_).N();

  if (activeVars_.has("streamfunction")) {
    applyNtimesNplus1SpectralScaling("vorticity", specFS_, N, fset["streamfunction"]);
  }
  if (activeVars_.has("velocity_potential")) {
    applyNtimesNplus1SpectralScaling("divergence", specFS_, N, fset["velocity_potential"]);
  }
  atlas::Field uvgp = allocateGaussUVField(outputFunctionSpace_,
                                           activeVars_, activeVariableSizes_);

  // transform to gaussian grid
  transFS_.invtrans_vordiv2wind(fset["vorticity"], fset["divergence"], uvgp);

  atlas::FieldSet uvfset = convertUVToFieldSet(uvgp);

  newFields.add(uvfset["eastward_wind"]);
  newFields.add(uvfset["northward_wind"]);

  fset = newFields;

  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}


// -----------------------------------------------------------------------------
void SPTOUV::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;

  // Create empty Model fieldset
  atlas::FieldSet newFields = atlas::FieldSet();

  // copy "passive variables"
  std::vector<std::string> fsetNames = fset.field_names();
  for (auto s : fset.field_names()) {
    if (!activeVars_.has(s)) {
      newFields.add(fset[s]);
    }
  }

  atlas::Field uvgp = convertUVToFieldSetAD(fset);

  atlas::FieldSet specfset = allocateSpectralVortDiv(specFS_, activeVars_, activeVariableSizes_);

  transFS_.invtrans_vordiv2wind_adj(uvgp, specfset["vorticity"], specfset["divergence"]);

  const int N = atlas::GaussianGrid(gaussGrid_).N();

  if (activeVars_.has("streamfunction")) {
    applyNtimesNplus1SpectralScaling("streamfunction", specFS_, N, specfset["vorticity"]);
  }
  if (activeVars_.has("velocity_potential")) {
    applyNtimesNplus1SpectralScaling("velocity_potential", specFS_, N, specfset["divergence"]);
  }

  newFields.add(specfset[0]);
  newFields.add(specfset[1]);

  fset = newFields;

  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------
void SPTOUV::calibrationInverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiplyAD starting" << std::endl;
  ABORT("SPTOUV::calibrationInverseMultiply: not implemented");
  oops::Log::trace() << classname() << "::inverseMultiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------
void SPTOUV::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace spectralb
}  // namespace saber
