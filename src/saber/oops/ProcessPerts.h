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
  typedef oops::IncrementMemberTemplateParameters<MODEL>     IncrementMemberTemplateParameters_;
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

  /// Where to read input ensemble: From states or perturbations
  oops::OptionalParameter<eckit::LocalConfiguration> ensemble{"ensemble", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensemblePert{"ensemble pert", this};

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
  typedef typename Increment_::WriteParameters_           IncrementWriteParameters_;
  typedef ProcessPertsParameters<MODEL>                   ProcessPertsParameters_;

 public:
// -----------------------------------------------------------------------------
  explicit ProcessPerts(const eckit::mpi::Comm & comm = eckit::mpi::comm()) :
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

    // Setup time
    const util::DateTime time = xx[0].validTime();

    eckit::LocalConfiguration filterCovarianceBlockConf(fullConfig, "saber filter blocks");
    std::unique_ptr<SaberParametricBlockChain> saberFilterBlocks;

    // List of output increments
    const auto & lowpassPerturbations = params.lowpassPerturbations.value();
    const IncrementWriteParameters_ & incrementsWriteParams =
      params.outputPerturbations;

    oops::Variables incVars = params.inputVariables;
    // Initialize outer variables
    const std::vector<std::size_t> vlevs = geom.variableSizes(incVars);
    for (std::size_t i = 0; i < vlevs.size() ; ++i) {
      incVars.addMetaData(incVars[i], "levels", vlevs[i]);
    }

    std::vector<oops::FieldSet3D> fsetEns;
    std::vector<oops::FieldSet3D> dualResFsetEns;
    eckit::LocalConfiguration covarConf;
    covarConf.set("iterative ensemble loading", false);
    covarConf.set("inverse test", false);
    covarConf.set("adjoint test", false);
    covarConf.set("square-root test", false);
    covarConf.set("covariance model", "SABER");
    covarConf.set("time covariance", "");

    // Initialize filter blockchain
    saberFilterBlocks = std::make_unique<SaberParametricBlockChain>(geom, geom,
      incVars, fsetXb, fsetFg, fsetEns, dualResFsetEns,
      covarConf, filterCovarianceBlockConf);

    // Yaml validation
    // TODO(Mayeul): Move this do an override of deserialize
    if (((params.ensemble.value() == boost::none) &&
        (params.ensemblePert.value() == boost::none)) ||
        ((params.ensemble.value() != boost::none) &&
        (params.ensemblePert.value() != boost::none)))
    {
      throw eckit::UserError(
       "Require either input states or input perturbations to be set in yaml",
       Here());
    }

    // Read input ensemble
    const bool iterativeEnsembleLoading = false;
    std::vector<oops::FieldSet3D> fsetEnsI;
    eckit::LocalConfiguration ensembleConf(fullConfig);
    readEnsemble<MODEL>(geom,
                        incVars,
                        xx[0], xx[0],
                        ensembleConf,
                        iterativeEnsembleLoading,
                        fsetEnsI);
    int nincrements = fsetEnsI.size();

    //  Loop over perturbations
    for (int jm = 0; jm < nincrements; ++jm) {
      //  Read ensemble member perturbation
      oops::FieldSet3D fsetI(fsetEnsI[jm]);
      Increment_ dxI(geom, incVars, time);
      dxI.zero();
      dxI.fromFieldSet(fsetI.fieldSet());

      //  Copy perturbation
      oops::FieldSet3D fset(fsetI);

      oops::Log::test() << "Norm of perturbation : member  " << jm+1
                        << ": " << dxI.norm() << std::endl;

      oops::FieldSet4D fset4dDxI(fsetI);
      oops::FieldSet4D fset4dDx(fset);

      // Apply filter blocks
      saberFilterBlocks->filter(fset4dDx);

      if (lowpassPerturbations != boost::none) {
        Increment_ dxLowPass(geom, incVars, time);
        dxLowPass.zero();
        dxLowPass.fromFieldSet(fset4dDx[0].fieldSet());

        auto lowpassPerturbationsUpdated(*lowpassPerturbations);
        lowpassPerturbationsUpdated.setMember(jm+1);
        dxLowPass.write(lowpassPerturbationsUpdated);
        oops::Log::test() << "Norm of low pass perturbation : member  " << jm+1
                          << ": " << dxLowPass.norm() << std::endl;
      }

      // High pass = full pert - low pass
      fset4dDx *= -1.0;
      fset4dDxI += fset4dDx;

      // Write high pass
      dxI.fromFieldSet(fset4dDxI[0].fieldSet());

      auto incrementsWriteParamsUpdated(incrementsWriteParams);
      incrementsWriteParamsUpdated.setMember(jm+1);
      dxI.write(incrementsWriteParamsUpdated);

      oops::Log::test() << "Norm of high pass perturbation : member  " << jm+1
                        << ": " << dxI.norm() << std::endl;
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
