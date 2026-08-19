/*
 * (C) Crown Copyright 2023-2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <algorithm>
#include <map>
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
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "saber/blocks/SaberParametricBlockChain.h"
#include "saber/oops/ErrorCovarianceParameters.h"
#include "saber/oops/Utilities.h"

namespace saber {

/// \brief 1) It first uses a vector of strings to successively dive into
///           the appropriate subconfiguration from a top level configuration
///           object.
///        2) It then takes this subconfiguration and reorders the keys
///           into a fixed order.
///        3) It then takes the values of each of the sorted keys and
///           concatenates them into a single string
///        4) All instances of the string "stringPattern" within the LocalConfiguration
///           conf are replaced by the concatenated string. The value "stringPattern"
///           is extracted from the configuration using the "patternNameKey".
void setConcatenatedString(const eckit::Configuration & fullConf,
                           const std::vector<std::string> & keyTags,
                           const std::string & patternNameKey,
                           eckit::LocalConfiguration & conf) {
  eckit::LocalConfiguration subconf(fullConf);
  for (const std::string& s : keyTags) {
    subconf = subconf.getSubConfiguration(s);
  }

  std::string vals("");
  std::vector<std::string> sortedKeys(subconf.keys());
  std::sort(sortedKeys.begin(), sortedKeys.end(),
            [](const std::string a, const std::string b) {return a > b; });

  for (const std::string& s : sortedKeys) {
    std::string val;
    subconf.get(s, val);
    vals.append(val);
  }

  if (conf.has(patternNameKey)) {
    const std::string stringPattern = conf.getString(patternNameKey);
    util::seekAndReplace(conf, stringPattern, vals);
  }
}

/// \brief Parameters for filtering a single perturbation
template <typename MODEL> class FilterParameters :
  public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(FilterParameters, oops::Parameters)

 public:
  /// Note that the parameters here are not actually used in the code
  /// They are here to express the intent of these variables.
  /// Later on in the code we use eckit::LocalConfiguration and check whether
  /// each key is in the configuration - if it is we use its value, otherwise
  /// it is set to "false".
  oops::Parameter<bool> residualFromFilter{
    "use residual from filter", false, this};

  // This is a vector of outer blocks defining the filter.
  oops::OptionalParameter<std::vector<SaberOuterBlockParametersWrapper>> filter{"filter", this};
};

// -----------------------------------------------------------------------------

/// \brief Write parameters for single filtered perturbation
template <typename MODEL> class OutputWriteParameters :
  public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(OutputWriteParameters, oops::Parameters)

 public:
  // This is a vector of outer blocks for diagnostic purposes.
  oops::OptionalParameter<std::vector<SaberOuterBlockParametersWrapper>> diagnosticOnlyBlock{
    "diagnostic only block", this};

  /// Write parameters using generic oops::util::writeFieldSet writer
  oops::OptionalParameter<eckit::LocalConfiguration>
    genericWrite{"generic write", this};

  /// Write parameters using model increment writer
  oops::OptionalParameter<eckit::LocalConfiguration>
    modelWrite{"model write", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class BandParameters :
  public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BandParameters, oops::Parameters)

 public:
  typedef FilterParameters<MODEL>                   FilterParameters_;
  typedef OutputWriteParameters<MODEL>              outputParameters_;

  oops::OptionalParameter<FilterParameters_> band{"band", this};
  oops::OptionalParameter<outputParameters_> output{"output", this};
};

// -----------------------------------------------------------------------------

/// \brief Top-level options taken by the ProcessPerts application.
template <typename MODEL> class ProcessPertsParameters :
  public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ProcessPertsParameters, oops::Parameters)

 public:
  typedef BandParameters<MODEL>                          BandParameters_;

  /// Geometry parameters.
  oops::RequiredParameter<eckit::LocalConfiguration> geometry{"geometry", this};

  /// Background parameters.
  oops::RequiredParameter<eckit::LocalConfiguration> background{"background", this};

  oops::RequiredParameter<oops::Variables> inputVariables{"input variables", this};

  oops::RequiredParameter<std::vector<BandParameters_>> bands{"bands", this};

  oops::Parameter<bool> recursivePertProcessing{"recursive perturbations processing", false, this};

  /// Where to read input ensemble: From states or perturbations
  oops::OptionalParameter<eckit::LocalConfiguration> ensemble{"ensemble", this};
  oops::OptionalParameter<eckit::LocalConfiguration> ensemblePert{"ensemble pert", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class ProcessPerts : public oops::Application {
  typedef oops::ModelSpaceCovarianceBase<MODEL>             CovarianceBase_;
  typedef oops::CovarianceFactory<MODEL>                    CovarianceFactory_;
  typedef oops::Geometry<MODEL>                             Geometry_;
  typedef oops::Increment<MODEL>                            Increment_;
  typedef oops::State<MODEL>                                State_;
  typedef oops::State4D<MODEL>                              State4D_;
  typedef ProcessPertsParameters<MODEL>                     ProcessPertsParameters_;

 public:
// -----------------------------------------------------------------------------
  explicit ProcessPerts(const eckit::mpi::Comm & comm = eckit::mpi::comm()) :
    Application(comm) {
    instantiateCovarFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~ProcessPerts() {}
// -----------------------------------------------------------------------------

  int execute(const eckit::Configuration & fullConfig) const override {
    // Define space and time communicators
    const eckit::mpi::Comm * commSpace = &this->getComm();
    const eckit::mpi::Comm * commTime = &oops::mpi::myself();

    // Replace patterns in full configuration and deserialize parameters
    eckit::LocalConfiguration fullConfigUpdated(fullConfig);
    const size_t ntasks = commSpace->size();
    size_t nthreads = 1;
#ifdef _OPENMP
    # pragma omp parallel
    {
      nthreads = omp_get_num_threads();
    }
#endif
    util::seekAndReplace(fullConfigUpdated, "_MPI_", std::to_string(ntasks));
    util::seekAndReplace(fullConfigUpdated, "_OMP_", std::to_string(nthreads));

    // Deserialize parameters
    ProcessPertsParameters_ params;
    params.deserialize(fullConfigUpdated);

    // Setup geometry
    const Geometry_ geom(params.geometry, *commSpace, *commTime);

    // Setup background
    const State4D_ xx(geom, params.background, *commTime);
    oops::FieldSet4D fsetXb(xx);
    oops::FieldSet4D fsetFg(xx);

    // Setup time
    const util::DateTime time = xx[0].validTime();

    oops::Variables incVars = params.inputVariables;

    // Initialize outer variables
    const std::vector<std::size_t> vlevs = geom.variableSizes(incVars);
    for (std::size_t i = 0; i < vlevs.size() ; ++i) {
      incVars[i].setLevels(vlevs[i]);
    }

    // Read input ensemble
    oops::FieldSets fsetEnsI = readEnsemble<MODEL>(geom,
                                                   incVars,
                                                   xx.times(), xx.commTime(), xx.commEns(),
                                                   fullConfigUpdated);
    int nincrements = fsetEnsI.ens_size();

    const std::size_t nbands = params.bands.value().size();
    const std::vector<eckit::LocalConfiguration> bandsConfs
      = fullConfigUpdated.getSubConfigurations("bands");
    const bool recursivePertProcessing = params.recursivePertProcessing.value();

    // need to create a vectors of saber block chains to use later
    std::map<std::size_t, std::vector<SaberOuterBlockParametersWrapper>> diagBlockConfs;
    std::map<std::size_t, std::vector<SaberOuterBlockParametersWrapper>> filterCovBlockConfs;
    std::map<std::size_t, eckit::LocalConfiguration> genericWriteConfs;
    std::map<std::size_t, eckit::LocalConfiguration> modelWriteConfs;
    std::vector<bool> calcComplement;

    std::size_t b(0);
    for (const auto & bandConf : bandsConfs) {
      if (bandConf.has("band")) {
        // Add filter for this band
        eckit::LocalConfiguration bConf = bandConf.getSubConfiguration("band");
        for (const auto & outerBlockConf : bConf.getSubConfigurations("filter")) {
          SaberOuterBlockParametersWrapper cmpOuterBlockParamsWrapper;
          cmpOuterBlockParamsWrapper.deserialize(outerBlockConf);
          filterCovBlockConfs[b].push_back(cmpOuterBlockParamsWrapper);
        }
        calcComplement.push_back(
          bConf.getBool("use residual from filter", false));
      } else {
        // Last band without filter: complement of the sum of all previous bands
        ASSERT(b == nbands-1);
      }

      if (bandConf.has("output")) {
        eckit::LocalConfiguration oConf = bandConf.getSubConfiguration("output");
        if (oConf.has("diagnostic only block")) {
          for (const auto & outerBlockConf : oConf.getSubConfigurations("diagnostic only block")) {
            SaberOuterBlockParametersWrapper cmpOuterBlockParamsWrapper;
            cmpOuterBlockParamsWrapper.deserialize(outerBlockConf);
            diagBlockConfs[b].push_back(cmpOuterBlockParamsWrapper);
          }
        }
        if (oConf.has("generic write")) {
          eckit::LocalConfiguration gConf = oConf.getSubConfiguration("generic write");
          genericWriteConfs[b] = gConf;
        }
        if (oConf.has("model write")) {
          eckit::LocalConfiguration mConf = oConf.getSubConfiguration("model write");
          modelWriteConfs[b] = mConf;
        }
      }
      b++;
    }

    std::vector<std::unique_ptr<SaberOuterBlockChain>> saberFilterBlocks;
    const ErrorCovarianceParametersBase paramsBase;
    for (const auto & [key, value] : filterCovBlockConfs) {
      saberFilterBlocks.push_back(
        std::make_unique<SaberOuterBlockChain>(geom,
                                               incVars,
                                               fsetXb,
                                               fsetFg,
                                               paramsBase.toConfiguration(),
                                               value));
    }

    std::vector<std::unique_ptr<SaberOuterBlockChain>> saberDiagnosticBlocks;
    for (const auto & [key, value] : diagBlockConfs) {
      saberDiagnosticBlocks.push_back(
        std::make_unique<SaberOuterBlockChain>(geom,
                                               incVars,
                                               fsetXb,
                                               fsetFg,
                                               paramsBase.toConfiguration(),
                                               value));
    }

    //  Loop over perturbations
    for (int jm = 0; jm < nincrements; ++jm) {
      // Initialize work perturbation xI from ensemble perturbation x0
      oops::FieldSet3D fsetI(fsetEnsI[jm]);

      oops::Log::test() << "Norm of perturbation: "
                        << "member " << jm+1
                        << ": " << fsetI.norm(fsetI.variables()) << std::endl;

      // Initialize sum of filtered perturbations
      oops::FieldSet3D fsetSum(fsetI.validTime(), fsetI.commGeom());
      fsetSum.allocateOnly(fsetI.fieldSet());
      fsetSum.zero();
      oops::FieldSet4D fset4dDxSum(fsetSum);

      for (std::size_t b = 0; b < nbands; ++b) {
        // Copy work perturbation x = xI
        oops::FieldSet3D fset(fsetI.validTime(), fsetI.commGeom());
        fset.deepCopy(fsetI.fieldSet());
        oops::FieldSet4D fset4dDx(fset);

        // Apply filter blocks
        if (auto it{filterCovBlockConfs.find(b)}; it != std::end(filterCovBlockConfs)) {
          // Get filter blocks index
          const std::size_t idx = std::distance(std::begin(filterCovBlockConfs), it);

          // Apply filter G on input x: x' = Gx
          saberFilterBlocks[idx]->applyOuterBlocks(fset4dDx);

          if (calcComplement[b]) {
            // Use filter complement: x' = (I-G)x
            fset4dDx[0] -= fsetI;
            fset4dDx[0] *= -1.0;
          }

          if (recursivePertProcessing) {
            // Recursive filters: xI = xI - x'
            fsetI -= fset4dDx[0];
          }

          // Increment sum with the latest x'
          fset4dDxSum += fset4dDx;
        } else {
          // Residual increment: x' = x0 - sum{previous x'}
          fset4dDx[0].zero();
          fset4dDx[0] += fsetEnsI[jm];
          fset4dDx[0] -= fset4dDxSum[0];
        }

        oops::Log::test() << "Norm of band perturbation: "
                          << "member " << jm+1 << ": band " << b+1
                          << ": " << fset4dDx[0].norm(fset4dDx[0].variables())
                          << std::endl;


        // Apply diagnostic blocks
        if (auto it{diagBlockConfs.find(b)}; it != std::end(diagBlockConfs)) {
          const std::size_t idx = std::distance(std::begin(diagBlockConfs), it);
          saberDiagnosticBlocks[idx]->applyOuterBlocksAD(fset4dDx);
          saberDiagnosticBlocks[idx]->applyOuterBlocks(fset4dDx);
        }

        if (auto it{genericWriteConfs.find(b)}; it != std::end(genericWriteConfs)) {
          eckit::LocalConfiguration gconf = it->second;
          util::setMember(gconf, jm+1);
          setConcatenatedString(fullConfigUpdated,
                                std::vector<std::string>{"geometry", "grid"},
                                "grid pattern",
                                gconf);
          util::writeFieldSet(geom.getComm(), gconf, fset4dDx[0].fieldSet());
        }

        if (auto it{modelWriteConfs.find(b)}; it != std::end(modelWriteConfs)) {
          eckit::LocalConfiguration mconf = it->second;

          // Should be on the model geometry!
          auto pert = Increment_(geom,
                                 fset4dDx[0].variables(),
                                 time);
          pert.zero();
          pert.fromFieldSet(fset4dDx[0].fieldSet());

          util::setMember(mconf, jm+1);
          pert.write(mconf);
        }
      }
    }

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::ProcessPerts<" + MODEL::name() + ">";
  }
};

}  // namespace saber
