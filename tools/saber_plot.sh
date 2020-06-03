#!/usr/bin/env bash
#----------------------------------------------------------------------
# Bash script: saber_plot
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Parameters
testdata=$1
test=$2
mpi=$3
omp=$4

# Available plots list
plot_list="cortrack"

# Activate conda
eval "$(conda shell.bash hook)"
conda activate pyn_env

# BUMP tests
if test "${test%%_*}" = "bump" ; then
   for file in `ls ${testdata}/${test}/test_${mpi}-${omp}_*.nc` ; do
      if test ! -L ${file}; then
         # Get suffix
         tmp=${file#${testdata}/${test}/test_${mpi}-${omp}_}
         suffix=${tmp%.nc}
         for plot in ${plot_list} ; do
            if printf %s\\n "${suffix}" | grep -qF "${plot}" ; then
               ./plot/saber_${suffix}.py ${testdata} ${test} test_${mpi}-${omp}_${suffix}
            fi
         done
      fi
   done
fi
