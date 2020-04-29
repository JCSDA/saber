/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 */

#ifndef TEST_SABER_BUMPINTERPOLATION_H_
#define TEST_SABER_BUMPINTERPOLATION_H_

#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/functionspace/PointCloud.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/option.h"
#include "atlas/util/Point.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "eckit/testing/Test.h"

#include "oops/generic/Interpolator.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"
#include "saber/interpolation/InterpolatorBump.h"
#include "oops/../test/TestEnvironment.h"

using atlas::array::make_view;
using atlas::option::halo;
using atlas::option::levels;
using atlas::option::name;

namespace saber {
namespace test {

// -----------------------------------------------------------------------------
/*! Test Fixture for interpolation interface */

class BumpInterpolationFixture : private boost::noncopyable {
 public:
  static const eckit::LocalConfiguration & test() {return *getInstance().test_;}

 private:
  static BumpInterpolationFixture & getInstance() {
  static BumpInterpolationFixture theBumpInterpolationFixture;
    return theBumpInterpolationFixture;
  }

  BumpInterpolationFixture() {
    test_.reset(new eckit::LocalConfiguration(::test::TestEnvironment::config(),
                                              "test_interpolation_interface"));
  }

  ~BumpInterpolationFixture() {}

  std::unique_ptr<const eckit::LocalConfiguration> test_;
};

// -----------------------------------------------------------------------------
/*! smooth function for testing interpolation */
double testfunc(double lon, double lat, std::size_t jlev = 1,
                                        std::size_t nlev = 1) {
  double deg2rad = M_PI / 180.0L;
  double zz = static_cast<double>(jlev+1)/static_cast<double>(nlev);
  return (sin(deg2rad*lon)*cos(deg2rad*lat) + zz)/1.5;
}

// -----------------------------------------------------------------------------
/*! Test C++ interface to Bump interpolation
 *
*/

void testBumpInterpolation() {
  typedef BumpInterpolationFixture Test_;

  // Skip this test if it's not specified in the yaml
  if (!Test_::test().has("bump_interpolation"))
    return;

  oops::Log::info() << "\nStarting test of oops bump  interpolation interface" << std::endl;

  // For testing use a regular Octahedral lat/lon grid
  atlas::Grid grid1("O64");

  // use 3 vertical levels for testing
  std::size_t nlev = 3;

  // Generate mesh and functionspace for input grid
  atlas::MeshGenerator meshgen("structured");
  atlas::Mesh mesh1 = meshgen.generate(grid1);
  atlas::functionspace::NodeColumns fs1(mesh1, levels(nlev) | halo(0));

  // Generate random point cloud as output grid
  eckit::LocalConfiguration config = Test_::test().getSubConfiguration("bump_interpolation");
  const std::size_t N = static_cast<size_t>(config.getInt("NRandom"));
  unsigned int seed = static_cast<unsigned int>(config.getInt("Seed"));

  // Distribute output longitudes over MPI tasks
  unsigned int check_no_obs = 0;
  if (config.has("check_no_obs")) {
    check_no_obs = static_cast<unsigned int>(config.getInt("check_no_obs"));
  }
  size_t ntasks = oops::mpi::comm().size();
  if (check_no_obs == 1) {
    ntasks = ntasks-1;
  }
  size_t myrank = oops::mpi::comm().rank();

  // number of output points on this mpi task
  size_t myN = 0;
  if (myrank < ntasks) {
    myN = N/ntasks;
    if (myrank < (N % ntasks)) myN += 1;
  }

  // longitude range assigned to this mpi task
  double mylon_min = myrank * 360.0/static_cast<double>(ntasks);
  double mylon_max = (myrank+1) * 360.0/static_cast<double>(ntasks);

  // task-dependent random seed
  unsigned int myseed = (myrank+1)*seed;

  util::UniformDistribution<double> lon(myN, mylon_min, mylon_max, myseed);
  util::UniformDistribution<double> lat(myN, -88.0, 88.0);
  lon.sort();

  atlas::PointXY point;
  std::vector<atlas::PointXY> gridpoints;
  for (std::size_t jj=0; jj < myN; ++jj) {
    point.assign(lon[jj], lat[jj]);
    gridpoints.push_back(point);
  }

  atlas::functionspace::PointCloud fs2(gridpoints);

  // construct bump interpolator object
  config.set("interpolator", "bump");
  config.set("nlevels", nlev);

  saber::InterpolatorBump interpolator(config, fs1, fs2);

  // Next - define the input fields
  atlas::Field field1 = fs1.createField<double>(name("testfield"));

  // now fill in the input field with a smooth test function
  auto lonlat = make_view<double, 2>(fs1.nodes().lonlat());
  auto infield = make_view<double, 2>(field1);

  for (size_t jnode = 0; jnode < static_cast<size_t>(fs1.nodes().size()); ++jnode) {
    for (size_t jlev = 0; jlev < static_cast<size_t>(fs1.levels()); ++jlev)
      infield(jnode, jlev) = testfunc(lonlat(jnode, 0), lonlat(jnode, 1), jlev, nlev);
  }

  atlas::FieldSet infields;
  infields.add(field1);

  // You can pass the bump interpolator an empty output field set and it will
  // allocate the output fields, compute them, and add them to the field set
  atlas::FieldSet outfields;

  // Apply the interpolation
  interpolator.apply(infields, outfields);

  oops::Log::info() << "\n----------------------------------------------------"
                    << "\nTesting bump interpolation"
                    << "\n----------------------------------------------------"
                    << std::endl;

  // get test tolerance
  const double tolerance = config.getDouble("tolerance");

  // check result from atlas interpolation
  auto lonlat2 = make_view<double, 2>(fs2.lonlat());
  auto outfield = make_view<double, 2>(outfields.field("testfield"));

  size_t jlev = 2;  // level for checking result

  for (size_t jnode = 0; jnode < static_cast<size_t>(fs2.size()); ++jnode) {
    EXPECT(oops::is_close(outfield(jnode, jlev),
           testfunc(lonlat2(jnode, 0), lonlat2(jnode, 1), jlev, nlev), tolerance));
  }

  oops::Log::info() << "\n----------------------------------------------------"
                    << "\nTesting bump interpolation adjoint"
                    << "\n----------------------------------------------------"
                    << std::endl;

  // define random field on output grid
  size_t nout = fs2.size()*nlev;
  util::UniformDistribution<double> x(nout, 0.0, 1.0);

  atlas::Field field2_adcheck = fs2.createField<double>(name("adcheck")|levels(nlev));
  auto adcheck2 = make_view<double, 2>(field2_adcheck);

  size_t idx = 0;
  for (size_t jnode = 0; jnode < static_cast<size_t>(fs2.size()); ++jnode) {
    for (size_t jlev = 0; jlev < static_cast<size_t>(nlev); ++jlev) {
      adcheck2(jnode, jlev) = x[idx];
      ++idx;
    }
  }

  atlas::FieldSet adcheck_grid2;
  adcheck_grid2.add(field2_adcheck);

  // Define empty FieldSet for the result (on grid1)
  atlas::FieldSet adcheck_grid1;

  // Apply interpolator adjoint
  interpolator.apply_ad(adcheck_grid2, adcheck_grid1);

  // Check adjoint

  // we need to implement a general dot product for atlas Fields.
  // For now, do it manually.

  double dot1 = 0;
  double mydot1 = 0;
  atlas::Field field1_adcheck = adcheck_grid1.field("adcheck");
  auto adcheck1 = make_view<double, 2>(field1_adcheck);
  atlas::mesh::IsGhostNode is_ghost(fs1.nodes());
  for (size_t jnode = 0; jnode < static_cast<size_t>(fs1.size()); ++jnode) {
    if (is_ghost(jnode) == 0) {
      for (size_t jlev = 0; jlev < nlev; ++jlev)
        mydot1 += infield(jnode, jlev)*adcheck1(jnode, jlev);
    }
  }
  oops::mpi::comm().allReduce(mydot1, dot1, eckit::mpi::sum());

  double dot2 = 0;
  double mydot2 = 0;
  for (size_t jnode = 0; jnode < static_cast<size_t>(fs2.size()); ++jnode) {
    for (size_t jlev = 0; jlev < nlev; ++jlev)
      mydot2 += outfield(jnode, jlev)*adcheck2(jnode, jlev);
  }
  oops::mpi::comm().allReduce(mydot2, dot2, eckit::mpi::sum());

  EXPECT(oops::is_close(dot1, dot2, tolerance));

  oops::Log::info() << "\nFinishing test of oops bump interpolation interface" << std::endl;
}

// -----------------------------------------------------------------------------

class BumpInterpolation : public oops::Test {
 public:
  BumpInterpolation() {}
  virtual ~BumpInterpolation() {}

 private:
  std::string testid() const {return "saber::test::BumpInterpolation";}

  void register_tests() const {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("saber/interpolation/testBumpInterpolation")
      { testBumpInterpolation(); });
  }
};

// -----------------------------------------------------------------------------
}  // namespace test
}  // namespace saber

#endif  // TEST_SABER_BUMPINTERPOLATION_H_
