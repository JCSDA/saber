#!/bin/sh
#----------------------------------------------------------------------
# Shell script: saber_compare
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Parameters
test=$1
mpi=$2
omp=$3

# Initialize exit status
status=0

# BUMP tests
if test "${test%%_*}" = "bump" ; then
   # NCCMP Parameters
   nthreads=4
   tolerance=1.e-4

   for file in `ls testdata/${test}/test_${mpi}-${omp}_*.nc` ; do
      if test ! -L ${file}; then
         # Get suffix
         tmp=${file#testdata/${test}/test_${mpi}-${omp}_}
         suffix=${tmp%.nc}

         if printf %s\\n "${suffix}" | grep -qF "distribution" ; then
            # Remove distribution
            rm -f ${file}
         else
            # Build reference file name, with a special case where the file is mpi-dependent
            fileref=testref/${test}/test_1-1_${suffix}.nc
            for special in "mom" "lct_cor" "nicas" "normality" "obs" "split" "vbal" ; do
               if printf %s\\n "${suffix}" | grep -qF "${special}" ; then
                  fileref=testref/${test}/test_${mpi}-1_${suffix}.nc
               fi
            done

            if [ -x "$(command -v nccmp)" ] ; then
               # Compare files with NCCMP
               echo "Command: nccmp -dfFmqS --threads=${nthreads} -T ${tolerance} ${file} ${fileref}"
               nccmp -dfFmqS --threads=${nthreads} -T ${tolerance} ${file} ${fileref}
               exit_code=$?
               if test "${exit_code}" != "0" ; then
                  echo "\e[31mTest failed (nccmp) checking: "${file#testdata/}"\e[0m"
                  status=1
                  exit ${status}
               fi
            else
               echo "\e[31mCannot find command: nccmp\e[0m"
               status=2
               exit ${status}
            fi
         fi
      fi
   done

   for file in `ls testoutput/${test}/test_${mpi}-${omp}.0000.out` ; do
      # Check for stars
      grep -q "*" ${file}
      exit_code=$?
      if test "${exit_code}" = "0" ; then
         echo "\e[31mTest failed (stars in output) checking: "${file#testoutput/}"\e[0m"
         status=3
         exit ${status}
      fi

      # Check for "Close listings" line
      grep -q "\-\-\- Close listings" ${file}
      exit_code=$?
      if test "${exit_code}" != "0" ; then
         echo "\e[31mTest failed (no listing closure) checking: "${file#testoutput/}"\e[0m"
         status=3
         exit ${status}
      fi
   done
fi

# Exit
exit ${status}
