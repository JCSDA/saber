/*!
 @brief Command line arguments parsing and call to the BUMP library
 @author Benjamin Menetrier
 @copyright This code is distributed under the CeCILL-C license
 @copyright Copyright 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT.
*/

#include <string.h>

#include "eckit/mpi/Comm.h"
#include "eckit/runtime/Main.h"

extern "C" {
  void bump_main_f90(const int &, const char *, const int &, const char *);
}

int main(int argc, char** argv) {
  if (argc != 3) return 1;
  eckit::Main::initialise(argc, argv);
  int n1 = strlen(argv[1]);
  int n2 = strlen(argv[2]);
  bump_main_f90(n1, argv[1], n2, argv[2]);
  eckit::mpi::finaliseAllComms();
  return 0;
}
