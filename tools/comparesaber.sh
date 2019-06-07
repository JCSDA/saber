#!/bin/sh
#----------------------------------------------------------------------
# Shell script: compare
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Parameters
test=$1
mpi=$2
omp=$3

# Initialize exit status
status=0

# BUMP tests
if test "${test%%_*}" = "bump" ; then
   # NCCMP Parameters
   nthreads=4
   tolerance=1.e-4

   for file in `ls testdata/${test}/test_${mpi}-${omp}_*.nc` ; do
      # Get suffix
      tmp=${file#testdata/${test}/test_${mpi}-${omp}_}
      suffix=${tmp%.nc}

      if printf %s\\n "${suffix}" | grep -qF "distribution" ; then
         # Remove distribution
         rm -f ${file}
      else
         # Build reference file name, with a special case where the file is mpi-dependent
         fileref=testref/${test}/test_1-1_${suffix}.nc
         for special in "mom" "nicas" "obs" "split" "vbal" ; do
            if printf %s\\n "${suffix}" | grep -qF "${special}" ; then
               fileref=testref/${test}/test_${mpi}-1_${suffix}.nc
            fi
         done

         if type "nccmp" > /dev/null ; then
            # Compare files with NCCMP
            echo "Command: nccmp -dfFmqS --threads=${nthreads} -T ${tolerance} ${file} ${fileref}"
            nccmp -dfFmqS --threads=${nthreads} -T ${tolerance} ${file} ${fileref}
            exit_code=$?
            if test "${exit_code}" != "0" ; then
               echo "\e[31mTest failed checking: "${file#testdata/}"\e[0m"
               status=1
            fi
         else
            echo "Cannot find command: nccmp"
         fi
      fi
   done

   for file in `ls testoutput/${test}/test_${mpi}-${omp}.0000.out` ; do
      # Check for stars
      grep -q "*" ${file}
      exit_code=$?
      if test "${exit_code}" = "0" ; then
         echo "\e[31mTest failed checking: "${file#testoutput/}"\e[0m"
         status=2
      fi

      # Check for "Close listings" line
      grep -q "\-\-\- Close listings" ${file}
            exit_code=$?
      if test "${exit_code}" != "0" ; then
         echo "\e[31mTest failed checking: "${file#testoutput/}"\e[0m"
         status=2
      fi
   done
fi

# Exit
exit ${status}
