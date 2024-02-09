/*
 * (C) Crown Copyright 2023-2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/FieldSets.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/State.h"
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

#include "saber/blocks/SaberParametricBlockChain.h"
#include "saber/oops/Utilities.h"

namespace saber {

std::string getGridId(const eckit::Configuration & fullConfig) {
  eckit::LocalConfiguration geomconf = fullConfig.getSubConfiguration("geometry");
  eckit::LocalConfiguration gridconf = geomconf.getSubConfiguration("grid");
  std::string gridtype;
  std::string N;
  std::string gridid;
  if (gridconf.has("name")) {
    gridconf.get("name", gridid);
  } else if (gridconf.has("type") && gridconf.has("N")) {
    gridconf.get("type", gridtype);
    gridconf.get("N", N);
    gridid = gridtype.append(N);
  } else {
    gridid = "_";
  }
  return gridid;
}

void setGridId(eckit::LocalConfiguration & conf, const std::string & gridid) {
  if (conf.has("grid id pattern")) {
    const std::string gridIdPattern = conf.getString("grid id pattern");
    util::seekAndReplace(conf, gridIdPattern, gridid);
  } else {
    conf.set("grid id pattern", gridid);
  }
}

// -----------------------------------------------------------------------------

/// \brief Write parameters for single filtered perturbation
template <typename MODEL> class ModelOrGenericWriteParameters :
  public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ModelOrGenericWriteParameters, oops::Parameters)

 public:
  typedef typename oops::Increment<MODEL>::WriteParameters_  IncrementWriteParameters_;

  /// Write parameters using generic oops::util::writeFieldSet writer
  oops::OptionalParameter<eckit::LocalConfiguration>
    genericWrite{"generic write", this};

  /// Write parameters using model increment writer
  oops::OptionalParameter<IncrementWriteParameters_>
    modelWrite{"model write", this};
};


// -----------------------------------------------------------------------------

/// \brief Top-level options taken by the ProcessPerts application.
template <typename MODEL> class ProcessPertsParameters :
  public oops::ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(ProcessPertsParameters, oops::ApplicationParameters)

 public:
  typedef oops::ModelSpaceCovarianceParametersWrapper<MODEL> CovarianceParameters_;

  typedef typename oops::Geometry<MODEL>::Parameters_     GeometryParameters_;
  typedef ErrorCovarianceParameters<MODEL>                ErrorCovarianceParameters_;
  typedef ModelOrGenericWriteParameters<MODEL>            ModelOrGenericWriteParameters_;
  /// Geometry parameters.
  oops::RequiredParameter<GeometryParameters_> geometry{"geometry", this};

  /// Background parameters.
  oops::RequiredParameter<eckit::LocalConfiguration> background{"background", this};

  oops::RequiredParameter<oops::Variables> inputVariables{"input variables", this};

  oops::RequiredParameter<std::vector<ErrorCovarianceParameters_>>
   saberFilterCovarianceParams{"saber filter blocks", this};

  /// Where to read input ensemble: From states or perturbations
  oops::OptionalParameter<eckit::LocalConfiguration> ensemble{"ensemble", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensemblePert{"ensemble pert", this};

  /// Optional residual filter
  oops::OptionalParameter<eckit::LocalConfiguration> residualFilter{"residual filter", this};

  /// Whether and where to write waveband perturbations
  oops::RequiredParameter<std::vector<ModelOrGenericWriteParameters_>>
    outputPerts{"output perturbations", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class ProcessPerts : public oops::Application {
  typedef oops::ModelSpaceCovarianceBase<MODEL>             CovarianceBase_;
  typedef oops::CovarianceFactory<MODEL>                    CovarianceFactory_;
  typedef oops::ModelSpaceCovarianceParametersBase<MODEL>   CovarianceParametersBase_;
  typedef oops::Geometry<MODEL>                             Geometry_;
  typedef oops::Increment<MODEL>                            Increment_;
  typedef oops::State<MODEL>                                State_;
  typedef oops::State4D<MODEL>                              State4D_;
  typedef typename oops::Increment<MODEL>::WriteParameters_ IncrementWriteParameters_;
  typedef ProcessPertsParameters<MODEL>                     ProcessPertsParameters_;

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

    std::vector<eckit::LocalConfiguration> filterCovarianceBlockConfs
      = fullConfig.getSubConfigurations("saber filter blocks");

    oops::Variables incVars = params.inputVariables;
    // Initialize outer variables
    const std::vector<std::size_t> vlevs = geom.variableSizes(incVars);
    for (std::size_t i = 0; i < vlevs.size() ; ++i) {
      incVars.addMetaData(incVars[i], "levels", vlevs[i]);
    }

    std::vector<util::DateTime> dates;
    std::vector<int> ensmems;
    oops::FieldSets fsetEns(dates, oops::mpi::myself(), ensmems, oops::mpi::myself());
    oops::FieldSets dualResFsetEns(dates, oops::mpi::myself(),
                                            ensmems, oops::mpi::myself());
    eckit::LocalConfiguration covarConf;
    covarConf.set("iterative ensemble loading", false);
    covarConf.set("inverse test", false);
    covarConf.set("adjoint test", false);
    covarConf.set("square-root test", false);
    covarConf.set("covariance model", "SABER");
    covarConf.set("time covariance", "");

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
    eckit::LocalConfiguration ensembleConf(fullConfig);
    eckit::LocalConfiguration outputEnsConf;
    oops::FieldSets fsetEnsI = readEnsemble<MODEL>(geom,
                                                   incVars,
                                                   xx, xx,
                                                   ensembleConf,
                                                   iterativeEnsembleLoading,
                                                   outputEnsConf);
    int nincrements = fsetEnsI.ens_size();

    std::vector<std::unique_ptr<SaberParametricBlockChain>> saberFilterBlocks;
    for (const auto & filterCovarianceBlockConf : filterCovarianceBlockConfs) {
      saberFilterBlocks.push_back(
        std::make_unique<SaberParametricBlockChain>(geom, geom,
                                                    incVars, fsetXb, fsetFg,
                                                    fsetEns, dualResFsetEns,
                                                    covarConf,
                                                    filterCovarianceBlockConf));
    }

    // Number of standard wavebands
    const std::size_t standardWavebands = filterCovarianceBlockConfs.size();

    // Optional additional waveband
    const auto & residualFilter = params.residualFilter.value();

    // Check that there are as many output perturbations as wavebands
    const std::size_t nOutputs = params.outputPerts.value().size();
    const std::size_t nWavebands = (residualFilter == boost::none) ?
                                    standardWavebands:
                                    standardWavebands + 1;
    if (nOutputs != nWavebands) {
      oops::Log::error() << "Error: Asked for " << nWavebands << " wavebands but "
                         << nOutputs << " output parameters were specified."
                         << std::endl;
      throw eckit::UserError("Number of wavebands and output perturbations inconsistent.",
                             Here());
    }

    //  Loop over perturbations
    for (int jm = 0; jm < nincrements; ++jm) {
      oops::FieldSet3D fsetI(fsetEnsI[jm]);
      oops::FieldSet4D fset4dDxI(fsetI);
      oops::Log::test() << "Norm of perturbation: "
                        << "member " << jm+1
                        << ": " << fsetI.norm(fsetI.variables()) << std::endl;

      oops::FieldSet3D fsetSum(fsetI.validTime(), fsetI.commGeom());
      fsetSum.allocateOnly(fsetI.fieldSet());
      fsetSum.zero();
      oops::FieldSet4D fset4dDxSum(fsetSum);

      for (std::size_t wb = 0; wb < standardWavebands; ++wb) {
        //  Copy perturbation
        oops::FieldSet3D fset(fsetI.validTime(), fsetI.commGeom());
        fset.deepCopy(fsetI.fieldSet());
        oops::FieldSet4D fset4dDx(fset);

        // Apply filter blocks
        saberFilterBlocks[wb]->filter(fset4dDx);

        fset4dDxSum += fset4dDx;

        // Optionally, write out perturbation
        auto localOutputPert = params.outputPerts.value()[wb];

        if (localOutputPert.genericWrite.value() != boost::none) {
          eckit::LocalConfiguration conf = localOutputPert.genericWrite.value().value();
          util::setMember(conf, jm+1);
          setGridId(conf, getGridId(fullConfig));
          util::writeFieldSet(geom.getComm(), conf, fset4dDx[0].fieldSet());
        }

        if (localOutputPert.modelWrite.value() != boost::none) {
          // Should be on the model geometry!
          auto pert = Increment_(geom,
                                 fset4dDx[0].variables(),
                                 time);
          pert.zero();
          pert.fromFieldSet(fset4dDx[0].fieldSet());

          IncrementWriteParameters_ writeParams = localOutputPert.modelWrite.value().value();
          writeParams.setMember(jm+1);
          pert.write(writeParams);
        }

        oops::Log::test() << "Norm of waveband perturbation: "
                          << "member " << jm+1 << ": waveband " << wb+1
                          << ": " << fset4dDx[0].norm(fset4dDx[0].variables()) << std::endl;
      }

      if (residualFilter != boost::none) {
        const std::size_t wb = standardWavebands;

        fset4dDxI[0] -= fset4dDxSum[0];

        // Optionally, write out perturbation
        auto localOutputPert = params.outputPerts.value()[wb];

        if (localOutputPert.genericWrite.value() != boost::none) {
          auto conf = localOutputPert.genericWrite.value().value();
          util::setMember(conf, jm+1);
          setGridId(conf, getGridId(fullConfig));
          util::writeFieldSet(geom.getComm(), conf, fset4dDxI[0].fieldSet());
        }

        if (localOutputPert.modelWrite.value() != boost::none) {
          // Should be on the model geometry!
          auto pert = Increment_(geom,
                                 fset4dDxI[0].variables(),
                                 time);
          pert.zero();
          pert.fromFieldSet(fset4dDxI[0].fieldSet());

          auto writeParams = localOutputPert.modelWrite.value().value();
          writeParams.setMember(jm+1);
          pert.write(writeParams);
        }

        oops::Log::test() << "Norm of waveband perturbation: "
                          << "member " << jm+1 << ": waveband " << wb+1
                          << ": " << fset4dDxI[0].norm(fset4dDxI[0].variables()) << std::endl;
      }
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
