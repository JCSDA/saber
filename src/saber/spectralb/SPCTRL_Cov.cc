/*
 * (C) Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "saber/spectralb/SPCTRL_Cov.h"

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"

#include "saber/oops/SaberCentralBlockBase.h"
#include "saber/oops/SaberCentralBlockParametersBase.h"
#include "saber/spectralb/spectralb.h"

namespace oops {
  class Variables;
}

namespace saber {
namespace spectralb {

// -----------------------------------------------------------------------------

static SaberCentralBlockMaker<SPCTRL_COV> makerSPCTRL_COV_("SPCTRL_COV");

// -----------------------------------------------------------------------------

SPCTRL_COV::SPCTRL_COV(const eckit::mpi::Comm & comm,
       const atlas::FunctionSpace & functionSpace,
       const atlas::FieldSet & extraFields,
       const std::vector<size_t> & variableSizes,
       const eckit::Configuration & conf,
       const atlas::FieldSet & xb,
       const atlas::FieldSet & fg,
       const std::vector<atlas::FieldSet> & fsetVec)
 : SaberCentralBlockBase(conf), spectralb_()
{
  oops::Log::trace() << classname() << "::SPCTRL_COV starting" << std::endl;

  // Deserialize configuration
  SPCTRL_COVParameters params;
  params.validateAndDeserialize(conf);

  // Initialize SpectralB
  spectralb_.reset(new SpectralB(functionSpace,
                                 variableSizes,
                                 *params.activeVars.value(),
                                 params.spectralbParams.value()));

  oops::Log::trace() << classname() << "::SPCTRL_COV done" << std::endl;
}

// -----------------------------------------------------------------------------

SPCTRL_COV::~SPCTRL_COV() {
  oops::Log::trace() << classname() << "::~SPCTRL_COV starting" << std::endl;
  util::Timer timer(classname(), "~SPCTRL_COV");
  oops::Log::trace() << classname() << "::~SPCTRL_COV done" << std::endl;
}

// -----------------------------------------------------------------------------

void SPCTRL_COV::randomize(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  ABORT("SPCTRL_COV::randomize: not implemented");
  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

void SPCTRL_COV::multiply(atlas::FieldSet & fset) const {
  oops::Log::trace() << classname() << "::multiply starting" << std::endl;
  spectralb_->multiply_InterpAndCov(fset);
  oops::Log::trace() << classname() << "::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void SPCTRL_COV::print(std::ostream & os) const {
  os << classname();
}

// -----------------------------------------------------------------------------

}  // namespace spectralb
}  // namespace saber
