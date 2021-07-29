/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef SABER_OOPS_ESTIMATEPARAMS_H_
#define SABER_OOPS_ESTIMATEPARAMS_H_

#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/Logger.h"

#include "saber/oops/instantiateCovarFactory.h"
#include "saber/oops/ParametersBUMP.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace saber {

template <typename MODEL> class EstimateParams : public oops::Application {
  typedef oops::Geometry<MODEL>                           Geometry_;
  typedef oops::Increment<MODEL>                          Increment_;
  typedef oops::State<MODEL>                              State_;
  typedef ParametersBUMP<MODEL>                           ParametersBUMP_;
  typedef oops::IncrementEnsemble<MODEL>                  Ensemble_;
  typedef std::shared_ptr<oops::IncrementEnsemble<MODEL>> EnsemblePtr_;

 public:
// -----------------------------------------------------------------------------
  static const std::string classname() {return "saber::EstimateParams";}
  explicit EstimateParams(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {
    instantiateCovarFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~EstimateParams() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
    util::Timer timer(classname(), "write");

    //  Setup resolution
    const eckit::LocalConfiguration resolConfig(fullConfig, "geometry");
    const Geometry_ resol(resolConfig, this->getComm());

    // Setup variables
    const oops::Variables vars(fullConfig, "input variables");

    // Setup background state
    const eckit::LocalConfiguration backgroundConfig(fullConfig, "background");
    State_ xx(resol, backgroundConfig);

    //  Setup time
    const util::DateTime time = xx.validTime();

    // Setup ensemble 1
    EnsemblePtr_ ens1 = NULL;
    if (fullConfig.has("ensemble")) {
      const eckit::LocalConfiguration ensembleConfig(fullConfig, "ensemble");
      ens1.reset(new Ensemble_(ensembleConfig, xx, xx, resol, vars));
      for (size_t ie = 0; ie < ens1->size(); ++ie) {
        (*ens1)[ie].toAtlas();
      }
    } else if (fullConfig.has("ensemble pert")) {
      const eckit::LocalConfiguration ensemblePertConfig(fullConfig, "ensemble pert");
      if (fullConfig.has("ensemble base")) {
        // Increment ensemble from difference of two state ensembles
        const eckit::LocalConfiguration ensembleBaseConfig(fullConfig, "ensemble base");
        ens1.reset(new Ensemble_(resol, vars, ensembleBaseConfig, ensemblePertConfig));
      } else {
        // Increment ensemble from increments on disk
        ens1.reset(new Ensemble_(resol, vars, ensemblePertConfig));
      }
      for (size_t ie = 0; ie < ens1->size(); ++ie) {
        (*ens1)[ie].toAtlas();
      }
    }

    // Setup ensemble 2
    EnsemblePtr_ ens2 = NULL;
    if (fullConfig.has("covariance")) {
      const eckit::LocalConfiguration covarConfig(fullConfig, "covariance");
      int ens2_ne = covarConfig.getInt("pseudoens_size");
      ens2.reset(new Ensemble_(resol, vars, time, ens2_ne));
      // One time-slot only
      std::unique_ptr<oops::ModelSpaceCovarianceBase<MODEL>>
        cov(oops::CovarianceFactory<MODEL>::create(covarConfig, resol, vars, xx, xx));
      for (int ie = 0; ie < ens2_ne; ++ie) {
        oops::Log::info() << "Generate pseudo ensemble member " << ie+1 << " / "
                          << ens2_ne << std::endl;

        // Compute a pseudo ensemble using randomization
        Increment_ incr(resol, vars, time);
        cov->randomize(incr);
        (*ens2)[ie] = incr;
        (*ens2)[ie].toAtlas();
      }
    }

    // Setup parameters
    ParametersBUMP_ param(resol, vars, fullConfig, ens1, ens2);

    // Write parameters
    param.write();

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "saber::EstimateParams<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace saber
#endif  // SABER_OOPS_ESTIMATEPARAMS_H_
