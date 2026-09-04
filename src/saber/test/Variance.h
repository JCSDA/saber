/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iomanip>
#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "atlas/field.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"

#include "oops/base/FieldSet3D.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment4D.h"
#include "oops/base/IncrementSet.h"
#include "oops/base/State4D.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/FloatCompare.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/IgnoreOtherParameters.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "test/TestEnvironment.h"

#include "saber/oops/ErrorCovariance.h"
#include "saber/oops/ErrorCovarianceParameters.h"

namespace saber {
namespace test {

// -----------------------------------------------------------------------------

/// Expected per-point variance for a variable: the global mean of the variance
/// over owned points and levels across all MPI tasks. A unit-diagonal
/// correlation operator gives 1.0; a StdDev outer block with standard deviation
/// sigma gives sigma^2.
class ExpectedVarianceNorm : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ExpectedVarianceNorm, oops::Parameters)
 public:
  oops::RequiredParameter<std::string> variable{"variable", this};
  oops::RequiredParameter<double> value{"value", this};
};

// -----------------------------------------------------------------------------

/// Optional Monte-Carlo cross-check of cov.variance() against the sample variance
/// of `samples` random draws from cov. Per-variable check on
///   ||cov.variance() - mc_variance|| / ||mc_variance|| < tolerance.
class MonteCarloCheck : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(MonteCarloCheck, oops::Parameters)
 public:
  oops::RequiredParameter<size_t> samples{"samples", this};
  oops::RequiredParameter<double> tolerance{"tolerance", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL>
class VarianceTestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(VarianceTestParameters, oops::Parameters)

 public:
  oops::RequiredParameter<eckit::LocalConfiguration> geometry{"geometry", this};
  oops::RequiredParameter<eckit::LocalConfiguration> background{"background", this};
  oops::RequiredParameter<ErrorCovarianceParameters> backgroundError{"background error", this};
  oops::OptionalParameter<std::vector<ExpectedVarianceNorm>> expected{
      "expected per-point variance", this};
  oops::Parameter<double> tolerance{"tolerance", 1.0e-12, this};
  oops::OptionalParameter<MonteCarloCheck> monteCarlo{"monte carlo", this};
  oops::IgnoreOtherParameters ignore{this};
};

// -----------------------------------------------------------------------------

template <typename MODEL>
class VarianceFixture : private boost::noncopyable {
 public:
  static const VarianceTestParameters<MODEL> & params() {return *getInstance().params_;}
  static void reset() {getInstance().params_.reset();}

 private:
  static VarianceFixture<MODEL> & getInstance() {
    static VarianceFixture<MODEL> inst;
    return inst;
  }

  VarianceFixture() {
    params_ = std::make_unique<VarianceTestParameters<MODEL>>();
    params_->validateAndDeserialize(::test::TestEnvironment::config());
  }

  ~VarianceFixture() {}

  std::unique_ptr<VarianceTestParameters<MODEL>> params_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
void testVariance() {
  using Geometry_     = oops::Geometry<MODEL>;
  using Increment4D_  = oops::Increment4D<MODEL>;
  using IncrementSet_ = oops::IncrementSet<MODEL>;
  using State4D_      = oops::State4D<MODEL>;

  const auto & params = VarianceFixture<MODEL>::params();

  const Geometry_ geom(params.geometry.value(), oops::mpi::world(), oops::mpi::myself());
  const State4D_ xx(geom, params.background.value(), oops::mpi::myself());
  const oops::Variables vars = xx.variables();

  saber::ErrorCovariance<MODEL> cov(geom, vars,
                                     params.backgroundError.value().toConfiguration(),
                                     xx, xx);

  oops::FieldSet3D variance = cov.variance();

  // Per-variable expected per-point variance, computed as the grid-size-
  // independent global mean
  //   mean(variance) = <variance, ones> / <ones, ones>.
  oops::FieldSet3D ones(xx[0].validTime(), geom.getComm());
  ones.init(variance[0].functionspace(), variance.variables(), 1.0);
  if (params.expected.value()) {
    for (const auto & expected : *params.expected.value()) {
      const oops::Variables thisvar(std::vector<std::string>{expected.variable.value()});
      const double sumVar = variance.dot_product_with(ones, thisvar);
      const double sumOne = ones.dot_product_with(ones, thisvar);
      const double computed = sumVar / sumOne;
      oops::Log::info() << std::setprecision(16)
                        << "Per-point variance for " << expected.variable.value()
                        << ": computed=" << computed
                        << " expected=" << expected.value.value() << std::endl;
      EXPECT(oops::is_close_relative(computed, expected.value.value(),
                                     params.tolerance.value()));
    }
  }

  // Monte-Carlo cross-check: compare cov.variance() to the sample variance of
  // `samples` random draws from cov, field-by-field, via a relative diff norm.
  if (params.monteCarlo.value()) {
    const size_t nSamples = params.monteCarlo.value()->samples.value();
    const double mcTol   = params.monteCarlo.value()->tolerance.value();
    ASSERT(nSamples > 1);

    // Build an IncrementSet with nSamples members on a single time slot, draw
    // each member from cov, then take the unbiased sample variance.
    std::vector<int> ensLabels(nSamples);
    for (size_t jm = 0; jm < nSamples; ++jm) ensLabels[jm] = static_cast<int>(jm);
    IncrementSet_ samples(geom, vars, std::vector<util::DateTime>{xx[0].validTime()},
                          oops::mpi::myself(), ensLabels, oops::mpi::myself());

    Increment4D_ scratch(geom, vars, std::vector<util::DateTime>{xx[0].validTime()},
                         oops::mpi::myself());
    for (size_t jm = 0; jm < nSamples; ++jm) {
      cov.randomize(scratch);
      samples(0, jm) = scratch[0];
    }

    const IncrementSet_ mcVar = samples.ens_var();
    oops::FieldSet3D mcVariance(xx[0].validTime(), geom.getComm());
    mcVariance.deepCopy(mcVar(0, 0).fieldSet());

    oops::Log::info() << "Monte-Carlo variance check with " << nSamples
                      << " samples and tolerance " << mcTol << std::endl;
    oops::Log::info() << "MC variance:" << std::endl << mcVariance << std::endl;
    oops::Log::info() << "Covariance variance:" << std::endl << variance << std::endl;

    // diff = cov.variance() - mcVariance, computed once; the norms below mask it
    // per variable.
    oops::FieldSet3D diff(xx[0].validTime(), geom.getComm());
    diff.deepCopy(variance);
    diff -= mcVariance;

    // Per-variable: ||cov.variance() - mcVariance|| / ||mcVariance|| < mcTol.
    for (const auto & var : vars) {
      const oops::Variables thisvar(std::vector<std::string>{var.name()});
      const double diffNorm = diff.norm(thisvar);
      const double mcNorm   = mcVariance.norm(thisvar);
      // When mcNorm is zero, fall back to the absolute difference norm.
      const double rel = (mcNorm > 0.0) ? diffNorm / mcNorm : diffNorm;
      oops::Log::info() << std::setprecision(16)
                        << "MC variance check for " << var.name()
                        << ": ||cov.variance() - mc|| = " << diffNorm
                        << " ; ||mc|| = " << mcNorm
                        << " ; relative = " << rel
                        << " ; tolerance = " << mcTol
                        << std::endl;
      EXPECT(rel < mcTol);
    }
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
class Variance : public oops::Test {
 public:
  using oops::Test::Test;
  ~Variance() override {VarianceFixture<MODEL>::reset();}

 private:
  std::string testid() const override {return "saber::test::Variance<" + MODEL::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test> & ts = eckit::testing::specification();
    ts.emplace_back(CASE("saber/Variance/testVariance") {testVariance<MODEL>();});
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test
}  // namespace saber
