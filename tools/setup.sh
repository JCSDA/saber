#!/bin/sh
#----------------------------------------------------------------------
# Shell script: setup
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Parameters
src=$1
dst=$2
test=$3
mpi=$4
omp=$5

# BUMP tests
if test "${test%%_*}" = "bump" ; then
   # Link to ensemble members
   for member in `ls ${src}/testdata/ens*.nc`; do
      member_base="$(basename -- $member)"
      ln -sf ${member} ${dst}/testdata/${test}/${member_base}
   done

   # Link to grid
   ln -sf ${src}/testdata/grid.nc ${dst}/testdata/${test}/grid.nc

   # Link wind fields
   for wind in `ls ${src}/testdata/wind*.nc`; do
      wind_base="$(basename -- $wind)"
      ln -sf ${wind} ${dst}/testdata/${test}/${wind_base}
   done

   # Set namelist
   sed -e s/_MPI_/${mpi}/g -e s/_OMP_/${omp}/g ${src}/testinput/${test}.yaml > ${dst}/testinput/${test}_${mpi}-${omp}.yaml
fi
