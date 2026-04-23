/*
 * (C) Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>

#include "eckit/mpi/Comm.h"

#include "oops/mpi/mpi.h"
#include "oops/runs/Run.h"

#include "saber/oops/ErrorCovarianceToolbox.h"
#include "saber/oops/instantiateCovarFactory.h"
#include "src/Traits.h"



int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  const eckit::mpi::Comm & initialDefaultComm = oops::mpi::world();

  // This application runs on one less than the comm world total of MPI ranks
  // for the ErrorCovarianceToolbox.
  // The final MPI rank is left alone.
  // This is done to simulate the case where we have an IO server
  // and where the executable does not see all the ranks.
  const size_t globalRank = initialDefaultComm.rank();
  const size_t globalSize = initialDefaultComm.size();

  const size_t myComponent = (globalRank >= globalSize -1 ? 1 : 0);

  const auto spaceCommName = ("default_space_" + std::to_string(myComponent));
  const auto & localSpaceComm = initialDefaultComm.split(myComponent, spaceCommName.c_str());
  eckit::mpi::setCommDefault(localSpaceComm.name().c_str());

  int stat = 0;
  if (myComponent == 0) {
    saber::instantiateCovarFactory<quench::Traits>();
    saber::ErrorCovarianceToolbox<quench::Traits> ect;
    stat = run.execute(ect, localSpaceComm);
  }
  return stat;
}
