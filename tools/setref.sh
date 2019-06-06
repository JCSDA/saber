#!/bin/sh
#----------------------------------------------------------------------
# Shell script: setref
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Parameters
testref=$1
testdata=$2
test=$3

# BUMP tests
if test "${test%%_*}" = "bump" ; then
   # Make directory
   cd ${testref}
   mkdir -p ${test}
   rm -f ${test}/*

   # Copy 1-1 files
   cp ${testdata}/${test}/test_1-1_*.nc ${test}

   # Copy 2-1 special files
   for special in "mom" "nicas" "obs" "split" "vbal" ; do
      if ls ${testdata}/${test}/test_2-1_${special}*.nc 1> /dev/null 2>&1; then
         cp -f ${testdata}/${test}/test_2-1_${special}*.nc ${test}
      fi
   done
fi
