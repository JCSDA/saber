/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/range/adaptors.hpp>

#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/oops/ReadInputFields.h"
#include "saber/oops/SaberBlockParametersBase.h"
#include "saber/oops/SaberCentralTBlock.h"
#include "saber/oops/SaberOuterTBlock.h"

namespace saber {

// -----------------------------------------------------------------------------

template <typename MODEL> class SaberBlockTestParameters
  : public oops::ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(SaberBlockTestParameters, oops::ApplicationParameters)

 public:
  typedef typename oops::Geometry<MODEL>::Parameters_ GeometryParameters_;
  typedef typename oops::State<MODEL>::Parameters_    StateParameters_;

  /// Geometry parameters
  oops::RequiredParameter<GeometryParameters_> geometry{"geometry", this};

  /// Tested variables
  oops::RequiredParameter<oops::Variables> variables{"variables", this};

  /// Background parameters
  oops::RequiredParameter<StateParameters_> background{"background", this};

  /// SABER blocks
  oops::OptionalParameter<SaberCentralTBlockParameters<MODEL>>
    saberCentralTBlockParams{"saber central block", this};
  oops::OptionalParameter<std::vector<SaberOuterTBlockParameters<MODEL>>>
    saberOuterTBlocksParams{"saber outer blocks", this};

  /// Adjoint test tolerance
  oops::Parameter<double> adjointTolerance{"adjoint test tolerance", 1.0e-12, this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class SaberBlockTest : public oops::Application {
  typedef oops::Geometry<MODEL>                                Geometry_;
  typedef oops::Increment<MODEL>                               Increment_;
  typedef typename boost::ptr_vector<SaberOuterTBlock<MODEL>>  SaberOuterTBlockVec_;
  typedef oops::State<MODEL>                                   State_;

 public:
  static const std::string classname() {return "saber::SaberBlockTest";}
  explicit SaberBlockTest(const eckit::mpi::Comm & comm = oops::mpi::world())
    : Application(comm) {}

  virtual ~SaberBlockTest() {}

  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
    util::Timer timer(classname(), "execute");

    // Deserialize parameters
    SaberBlockTestParameters<MODEL> params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);

    // Setup geometry
    const Geometry_ geom(params.geometry, this->getComm());

    // Setup variables
    const oops::Variables incVars(params.variables);

    // Setup background state
    const State_ xx(geom, params.background);

    // SABER blocks
    std::unique_ptr<SaberCentralTBlock<MODEL>> saberCentralTBlock_;
    SaberOuterTBlockVec_ saberOuterTBlocks_;

    // Local copy of background and first guess
    State_ xbLocal(xx);
    State_ fgLocal(xx);

    // Initial outer geometry and variables
    std::vector<std::reference_wrapper<const oops::GeometryData>> outerGeometryData_;
    outerGeometryData_.push_back(geom.generic());
    oops::Variables outerVars(incVars);

    // Build outer blocks successively
    const boost::optional<std::vector<SaberOuterTBlockParameters<MODEL>>> &saberOuterTBlocksParams =
      params.saberOuterTBlocksParams.value();
    if (saberOuterTBlocksParams != boost::none) {
      // Loop in reverse order
      for (const SaberOuterTBlockParameters<MODEL> & saberOuterTBlockParams :
        boost::adaptors::reverse(*saberOuterTBlocksParams)) {

        // Create outer templated block
        saberOuterTBlocks_.push_back(new SaberOuterTBlock<MODEL>(geom,
                                     outerGeometryData_.back().get(),
                                     outerVars,
                                     saberOuterTBlockParams,
                                     xbLocal,
                                     fgLocal));

        // Update outer geometry and variables for the next block
        outerGeometryData_.push_back(saberOuterTBlocks_.back().innerGeometryData());
        outerVars = saberOuterTBlocks_.back().innerVars();
      }
    }

    // Get central templated block parameters
    const boost::optional<SaberCentralTBlockParameters<MODEL>> &saberCentralTBlockParams =
      params.saberCentralTBlockParams.value();

    if (saberCentralTBlockParams != boost::none) {
      // Create central templated block
      saberCentralTBlock_.reset(new SaberCentralTBlock<MODEL>(geom,
                                outerGeometryData_.back().get(),
                                outerVars,
                                *saberCentralTBlockParams,
                                xbLocal,
                                fgLocal,
                                params.adjointTolerance.value()));
    }

    return 0;
  }

 private:
  std::string appname() const override {
    return "saber::SaberBlockTest<" + MODEL::name() + ">";
  }
};

// -----------------------------------------------------------------------------

}  // namespace saber
