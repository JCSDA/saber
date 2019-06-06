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
   # Link to ensemble 1
   for timeslot in "00" "06" ; do
      for member in $(seq -f "%04g" 1 50) ; do
         ln -sf ${src}/testdata/ens1_${timeslot}_${member}.nc ${dst}/testdata/${test}/ens1_${timeslot}_${member}.nc
      done
   done

   # Link to ensemble 2
   for timeslot in "00" "06" ; do
      for member in $(seq -f "%04g" 1 10) ; do
         ln -sf ${src}/testdata/ens2_00_${member}.nc ${dst}/testdata/${test}/ens2_${timeslot}_${member}.nc
      done
   done

   # Link to grid
   ln -sf ${src}/testdata/grid.nc ${dst}/testdata/${test}/grid.nc

   # Set namelist
   sed -e s/_MPI_/${mpi}/g -e s/_OMP_/${omp}/g ${src}/testinput/${test}.nam > ${dst}/testinput/${test}_${mpi}-${omp}.nam
fi
