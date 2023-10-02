/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <omp.h>

#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Geometry.h"

#include "oops/base/Increment.h"
#include "oops/base/Increment4D.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/State.h"
#include "oops/base/State4D.h"
#include "oops/base/StateEnsemble.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/DateTime.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/blocks/SaberBlockChainBase.h"
#include "saber/blocks/SaberBlockParametersBase.h"
#include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/blocks/SaberOuterBlockChain.h"
#include "saber/blocks/SaberParametricBlockChain.h"

#include "saber/oops/Utilities.h"

namespace saber {

// -----------------------------------------------------------------------------

/// \brief Top-level options taken by the ProcessPerts application.
template <typename MODEL> class ProcessPertsParameters :
  public oops::ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(ProcessPertsParameters, oops::ApplicationParameters)

 public:
  typedef oops::ModelSpaceCovarianceParametersWrapper<MODEL> CovarianceParameters_;

  typedef typename oops::Geometry<MODEL>::Parameters_        GeometryParameters_;
  typedef oops::State<MODEL>                                 State_;
  typedef oops::StateEnsembleParameters<MODEL>               StateEnsembleParameters_;
  typedef typename oops::Increment<MODEL>::ReadParameters_   IncrementReadParameters_;
  typedef typename oops::Increment<MODEL>::WriteParameters_  IncrementWriterParameters_;
  typedef ErrorCovarianceParameters<MODEL>                   ErrorCovarianceParameters_;

  /// Geometry parameters.
  oops::RequiredParameter<GeometryParameters_> geometry{"geometry", this};

  /// Background parameters.
  oops::RequiredParameter<eckit::LocalConfiguration> background{"background", this};

  oops::RequiredParameter<util::DateTime> date{"date", this};
  oops::RequiredParameter<oops::Variables> inputVariables{"input variables", this};

  oops::RequiredParameter<ErrorCovarianceParameters_>
    saberFilterCovarianceParams{"saber filter blocks", this};

  /// Where to read optional input perturbations
  oops::OptionalParameter<std::vector<IncrementReadParameters_>>
    inputPerturbations{"input perturbations", this};

  /// Where to read optional input states
  oops::OptionalParameter<StateEnsembleParameters_>
    ensembleParams{"ensemble", this};

  /// Where to write low-pass filtered perturbations
  oops::OptionalParameter<IncrementWriterParameters_>
    lowpassPerturbations{"low pass perturbations", this};

  /// Where to write high-pass filtered perturbations
  oops::RequiredParameter<IncrementWriterParameters_>
    outputPerturbations{"output perturbations", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class ProcessPerts : public oops::Application {
  typedef oops::ModelSpaceCovarianceBase<MODEL>           CovarianceBase_;
  typedef oops::CovarianceFactory<MODEL>                  CovarianceFactory_;
  typedef oops::ModelSpaceCovarianceParametersBase<MODEL> CovarianceParametersBase_;
  typedef oops::Geometry<MODEL>                           Geometry_;
  typedef oops::Increment<MODEL>                          Increment_;
  typedef oops::Increment4D<MODEL>                        Increment4D_;
  typedef oops::State<MODEL>                              State_;
  typedef oops::State4D<MODEL>                            State4D_;
  typedef oops::StateEnsemble<MODEL>                      StateEnsemble_;
  typedef oops::StateEnsembleParameters<MODEL>            StateEnsembleParameters_;
  typedef typename Increment_::ReadParameters_            IncrementReadParameters_;
  typedef typename Increment_::WriteParameters_           IncrementWriteParameters_;
  typedef ProcessPertsParameters<MODEL>                   ProcessPertsParameters_;

 public:
// -----------------------------------------------------------------------------
  explicit ProcessPerts(const eckit::mpi::Comm & comm = oops::mpi::world()) :
    Application(comm) {
    oops::instantiateCovarFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~ProcessPerts() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
    // Deserialize parameters
    ProcessPertsParameters_ params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);

    // Define number of subwindows
    const eckit::LocalConfiguration backgroundConfig(fullConfig, "background");

    // Define space and time communicators
    const eckit::mpi::Comm * commSpace = &this->getComm();
    const eckit::mpi::Comm * commTime = &oops::mpi::myself();

    // Setup geometry
    const Geometry_ geom(params.geometry, *commSpace, *commTime);

    // Setup background
    const State4D_ xx(geom, params.background, *commTime);
    oops::FieldSet4D fsetXb(xx);
    oops::FieldSet4D fsetFg(xx);

    // Setup variables
    const oops::Variables statevars = xx.variables();

    // Setup time
    const util::DateTime time = xx[0].validTime();

    eckit::LocalConfiguration filterCovarianceBlockConf(fullConfig, "saber filter blocks");
    std::unique_ptr<SaberParametricBlockChain> saberFilterBlocks;

    //  List of input and output increments
    const auto & incrementsReadParams = params.inputPerturbations.value();
    const auto & ensembleParams = params.ensembleParams.value();

    // List of output increments
    const auto & lowpassPerturbations = params.lowpassPerturbations.value();
    const IncrementWriteParameters_ & incrementsWriteParams =
      params.outputPerturbations;

    const oops::Variables incvars = params.inputVariables;

    std::vector<atlas::FieldSet> fsetEns;
    std::vector<atlas::FieldSet> dualResolutionFsetEns;
    eckit::LocalConfiguration covarConf;
    covarConf.set("iterative ensemble loading", false);
    covarConf.set("inverse test", false);
    covarConf.set("adjoint test", false);
    covarConf.set("covariance model", "SABER");
    covarConf.set("time covariance", "");

    // Initialize filter blockchain
    saberFilterBlocks = std::make_unique<SaberParametricBlockChain>(geom, geom,
      incvars, fsetXb, fsetFg, fsetEns, dualResolutionFsetEns,
      covarConf, filterCovarianceBlockConf);

    int nincrements(0);
    if (incrementsReadParams != boost::none) {
      nincrements = (*incrementsReadParams).size();
    }

    if (((incrementsReadParams == boost::none) &&
        (ensembleParams == boost::none)) ||
        ((incrementsReadParams != boost::none) &&
        (ensembleParams != boost::none)))
    {
      throw eckit::UserError(
       "Require either input states or input perturbations to be set in yaml",
       Here());
    }

    // create ensemble mean states if states are read in
    State_ meanState(geom, incvars, time);
    meanState.zero();
    if (ensembleParams != boost::none) {
      nincrements = (*ensembleParams).size();
    }
    StateEnsemble_ stateEnsemble(meanState, nincrements);
    if (ensembleParams != boost::none) {
      StateEnsemble_ stateEnsembleTmp(geom, *ensembleParams);
      for (int jm = 0; jm < nincrements; ++jm) {
        util::copyFieldSet(stateEnsembleTmp[jm].fieldSet(),
                           stateEnsemble[jm].fieldSet());
      }
    }
    util::copyFieldSet(stateEnsemble.mean().fieldSet(),
                       meanState.fieldSet());

    //  Loop over perturbations
    for (int jm = 0; jm < nincrements; ++jm) {
      //  Read ensemble member perturbations
      Increment_ dxI(geom, incvars, time);

      if (incrementsReadParams != boost::none) {
        dxI.read((*incrementsReadParams)[jm]);
      }
      if (ensembleParams != boost::none) {
        dxI.diff(stateEnsemble[jm], meanState);
      }

      //  Copy pert to give #Copy
      Increment_ dx(dxI, true);
      oops::FieldSet4D fset4dDxI(oops::FieldSet3D{dxI.fieldSet(), time, geom.getComm()});
      oops::FieldSet4D fset4dDx(oops::FieldSet3D{dx.fieldSet(), time, geom.getComm()});

      // Apply filter blocks
      saberFilterBlocks->filter(fset4dDx);

      if (lowpassPerturbations != boost::none) {
        Increment_ dxLowPass(geom, incvars, time);
        dxLowPass.zero();
        dxLowPass.fromFieldSet(fset4dDx[0].fieldSet());

        auto lowpassPerturbationsUpdated(*lowpassPerturbations);
        lowpassPerturbationsUpdated.setMember(jm+1);
        dxLowPass.write(lowpassPerturbationsUpdated);
      }

      // High pass = full pert - low pass
      fset4dDx *= -1.0;
      fset4dDxI += fset4dDx;

      // Write #Filtered.
      Increment_ dxO(geom, incvars, time);
      dxO.zero();
      dxO.fromFieldSet(fset4dDxI[0].fieldSet());

      auto incrementsWriteParamsUpdated(incrementsWriteParams);
      incrementsWriteParamsUpdated.setMember(jm+1);
      dxO.write(incrementsWriteParamsUpdated);
    }

    return 0;
  }
// -----------------------------------------------------------------------------
  void outputSchema(const std::string & outputPath) const override {
    ProcessPertsParameters_ params;
    params.outputSchema(outputPath);
  }
// -----------------------------------------------------------------------------
  void validateConfig(const eckit::Configuration & fullConfig) const override {
    ProcessPertsParameters_ params;
    params.validate(fullConfig);
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::ProcessPerts<" + MODEL::name() + ">";
  }
};

}  // namespace saber
