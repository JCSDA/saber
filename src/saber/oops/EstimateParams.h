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
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/assimilation/Increment4D.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/parallel/mpi/mpi.h"
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
  typedef oops::Increment4D<MODEL>                        Increment4D_;
  typedef oops::State<MODEL>                              State_;
  typedef oops::State4D<MODEL>                            State4D_;
  typedef ParametersBUMP<MODEL>                           Parameters_;
  typedef oops::IncrementEnsemble<MODEL>                  Ensemble_;
  typedef std::shared_ptr<oops::IncrementEnsemble<MODEL>> EnsemblePtr_;

 public:
// -----------------------------------------------------------------------------
  static const std::string classname() {return "saber::EstimateParams";}
  explicit EstimateParams(const eckit::mpi::Comm & comm = oops::mpi::comm()) : Application(comm) {
    instantiateCovarFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~EstimateParams() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
    util::Timer timer(classname(), "write");

    //  Setup resolution
    const eckit::LocalConfiguration resolConfig(fullConfig, "resolution");
    const Geometry_ resol(resolConfig, this->getComm());

    // Setup variables
    const oops::Variables vars(fullConfig);

    // Setup background state
    const eckit::LocalConfiguration backgroundConfig(fullConfig, "background");
    State4D_ xx(resol, vars, backgroundConfig);

    //  Setup timeslots
    std::vector<util::DateTime> timeslots = xx.validTimes();
    oops::Log::info() << "Number of ensemble time-slots:" << timeslots.size() << std::endl;

    // Setup ensemble
    EnsemblePtr_ ens = NULL;
    if (fullConfig.has("ensemble")) {
      const eckit::LocalConfiguration ensembleConfig(fullConfig, "ensemble");
      ens.reset(new Ensemble_(ensembleConfig, xx, xx, resol));
    }

    // Setup pseudo ensemble
    EnsemblePtr_ pseudo_ens = NULL;
    if (fullConfig.has("covariance")) {
      const eckit::LocalConfiguration covarConfig(fullConfig, "covariance");
      int ens2_ne = covarConfig.getInt("pseudoens_size");
      pseudo_ens.reset(new Ensemble_(resol, vars, timeslots, ens2_ne));
      if (timeslots.size() == 1) {
      // One time-slot only
        std::unique_ptr<oops::ModelSpaceCovarianceBase<MODEL>>
          cov(oops::CovarianceFactory<MODEL>::create(covarConfig, resol, vars, xx[0], xx[0]));
        for (int ie = 0; ie < ens2_ne; ++ie) {
          oops::Log::info() << "Generate pseudo ensemble member " << ie+1 << " / "
                      << ens2_ne << std::endl;

          // Compute a pseudo ensemble using randomization
          Increment4D_ incr(resol, vars, timeslots);
          cov->randomize(incr[incr.first()]);
          (*pseudo_ens)[ie] = incr;
        }
      } else {
        // Multiple time-slots
        std::unique_ptr<oops::ModelSpaceCovariance4DBase<MODEL>>
          cov(oops::Covariance4DFactory<MODEL>::create(covarConfig, resol, vars, xx, xx));
        for (int ie = 0; ie < ens2_ne; ++ie) {
          oops::Log::info() << "Generate pseudo ensemble member " << ie+1 << " / "
                      << ens2_ne << std::endl;

          // Compute a pseudo ensemble using randomization
          Increment4D_ incr(resol, vars, timeslots);
          cov->randomize(incr);
          (*pseudo_ens)[ie] = incr;
        }
      }
    }

    // Setup parameters
    Parameters_ param(resol, vars, timeslots, fullConfig, ens, pseudo_ens);

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
