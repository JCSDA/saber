#!/bin/sh
#----------------------------------------------------------------------
# Shell script: saber_compare_oops
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Parameters
test=$1

# Grep
grep 'Test     : ' testoutput/${test}/test.log.out > testdata/${test}/test.log.out

# Diff
diff -s testref/${test}/test.log.out testdata/${test}/test.log.out
