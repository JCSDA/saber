/*!
 @brief Command line arguments parsing and call to the BUMP library
 @author Benjamin Menetrier
 @copyright This code is distributed under the CeCILL-C license
 @copyright Copyright 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT.
*/

#include <string.h>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/config/YAMLConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "eckit/runtime/Main.h"

extern "C" {
  void bump_main_f90(const eckit::Configuration &, const int &, const char *);
}

int main(int argc, char** argv) {
  if (argc != 3) return 1;
  eckit::Main::initialise(argc, argv);

  // Get configuration file
  eckit::PathName configfile = argv[1];

  // Read configuration
  const eckit::YAMLConfiguration config(configfile);

  // Print configuration
  std::cout << "Configuration input file is: " << configfile << std::endl;
  std::cout << "Full configuration is:"  << config << std::endl;

  // Call Fortran
  int n2 = strlen(argv[2]);
  bump_main_f90(config, n2, argv[2]);

  // Finalise communicators
  eckit::mpi::finaliseAllComms();
  return 0;
}
