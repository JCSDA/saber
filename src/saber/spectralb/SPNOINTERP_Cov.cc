/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/spectralb/SPNOINTERP_Cov.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"

#include "saber/oops/SaberBlockBase.h"
#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/spectralb/spectralbnointerp.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

static SaberBlockMaker<SPNOINTERP_COV>  makerSPNOINTERP_COV_("SPNOINTERP_COV");

// -----------------------------------------------------------------------------

SPNOINTERP_COV::SPNOINTERP_COV(const eckit::mpi::Comm & comm,
                                      const atlas::FunctionSpace & functionSpace,
                                      const atlas::FieldSet & extraFields,
                                      const std::vector<size_t> & variableSizes,
                                      const Parameters_ & params,
                                      const atlas::FieldSet & xb,
                                      const atlas::FieldSet & fg,
                                      const std::vector<atlas::FieldSet> & fsetVec)
  : SaberBlockBase(params), spectralb_()
{
  oops::Log::trace() << classname() << "::SPNOINTERP_COV starting" << std::endl;

  // Setup and check input/ouput variables
  const oops::Variables inputVars = params.inputVars.value();
  const oops::Variables outputVars = params.outputVars.value();
  ASSERT(inputVars == outputVars);

  // Active variables
  const boost::optional<oops::Variables> &activeVarsPtr = params.activeVars.value();
  oops::Variables activeVars;
  if (activeVarsPtr != boost::none) {
    activeVars += *activeVarsPtr;
    ASSERT(activeVars <= inputVars);
  } else {
    activeVars += inputVars;
  }

  // Initialize SpectralBNoInterp
  spectralb_.reset(new SpectralBNoInterp(variableSizes,
                                         activeVars,
                                         params.spectralbParams.value()));

  oops::Log::trace() << classname() << "::SPNOINTERP_COV done" << std::endl;
}

// -----------------------------------------------------------------------------

SPNOINTERP_COV::~SPNOINTERP_COV() {
  oops::Log::trace() << classname() << "::~SPNOINTERP_COV starting" << std::endl;
  util::Timer timer(classname(), "~SPNOINTERP_COV");
  oops::Log::trace() << classname() << "::~SPNOINTERP_COV done" << std::endl;
}

// -----------------------------------------------------------------------------

void SPNOINTERP_COV::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  ABORT("SPNOINTERP_COV::randomize: not implemented");
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void SPNOINTERP_COV::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  spectralb_->multiply(fset);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void SPNOINTERP_COV::inverseMultiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiply starting" << std::endl;
  ABORT("SPNOINTERP_COV::inverseMultiply: not implemented");
  oops::Log::trace() << classname() << "::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void SPNOINTERP_COV::multiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiplyAD starting" << std::endl;
  ABORT("SPNOINTERP_COV::multiplyAD: not implemented");
  oops::Log::trace() << classname() << "::multiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void SPNOINTERP_COV::inverseMultiplyAD(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::inverseMultiplyAD starting" << std::endl;
  ABORT("SPNOINTERP_COV::inverseMultiplyAD: not implemented");
  oops::Log::trace() << classname() << "::inverseMultiplyAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void SPNOINTERP_COV::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace spectralb
}  // namespace saber
