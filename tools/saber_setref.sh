#!/bin/sh
#----------------------------------------------------------------------
# Shell script: saber_setref
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Parameters
testref=$1
testdata=$2
test=$3

# Make directory
mkdir -p ${testdata}/${test}

# BUMP tests
if test "${test%%_*}" = "bump"; then
   # Copy 1-1 files
   for file in `ls ${testdata}/${test}/test_1-1_*.nc`; do
      if test ! -L ${file}; then
         echo "$(basename -- ${file})"
         cp -f ${file} ${testref}
      fi
   done

   # Copy 2-1 and 4-1 special files
   for special in "mom" "lct_cor" "nicas" "obs" "split" "vbal"; do
      if ls ${testdata}/${test}/test_2-1_${special}*.nc 1> /dev/null 2>&1; then
         for file in `ls ${testdata}/${test}/test_2-1_${special}*.nc`; do
            if test ! -L ${file}; then
               echo "$(basename -- $file)"
               cp -f ${file} ${testref}
            fi
         done
      fi
      if ls ${testdata}/${test}/test_4-1_${special}*.nc 1> /dev/null 2>&1; then
         for file in `ls ${testdata}/${test}/test_4-1_${special}*.nc`; do
            if test ! -L ${file}; then
               echo "$(basename -- $file)"
               cp -f ${file} ${testref}
            fi
         done
      fi
   done
fi
