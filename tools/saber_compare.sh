#!/usr/bin/env bash
#----------------------------------------------------------------------
# Bash script: saber_compare
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Parameters
test=$1

# Special suffixes lists
mpi_dependent="mom lct_cor nicas_grids nicas_local normality sampling_local sampling_grids vbal_cov_local vbal_local"
not_compared=""

# Initialize exit status
status=0

# BUMP tests
if test "${test%%_*}" = "bump" ; then
   # NCCMP Parameters
   nthreads=4
   tolerance=1.e-4

   # Get comparison type
   compare_type=''
   if test "$#" = 3 ; then
      # Arguments 2 and 3 are integer => normal
      if [[ $2 =~ ^-?[0-9]+$ ]] && [[ $3 =~ ^-?[0-9]+$ ]] ; then
         compare_type="normal"
      fi
      if [[ ! $2 =~ ^-?[0-9]+$ ]] && [[ ! $4 =~ ^-?[0-9]+$ ]] && [[ ! $3 =~ ^-?[0-9]+$ ]] ; then
         compare_type="specific"
      fi
   fi

   # Check comparison type
   if test "${compare_type}" = '' ; then
      echo -e "\e[31mWrong number/type of arguments\e[0m"
      exit
   fi

   if test "${compare_type}" = "normal" ; then
      # Normal tests
      mpi=$2
      omp=$3

      # Check number of files
      nfiles=`ls testdata/${test}/test_${mpi}-${omp}_*.nc 2>/dev/null | wc -l`
      if test "${nfiles}" = "0"; then
         echo -e "\e[31mNo NetCDF file to check\e[0m"
         status=1
      fi

      for file in `ls testdata/${test}/test_${mpi}-${omp}_*.nc 2>/dev/null` ; do
         if test ! -h ${file}; then
            # Get suffix
            tmp=${file#testdata/${test}/test_${mpi}-${omp}_}
            suffix=${tmp%.nc}

            if printf %s\\n "${suffix}" | grep -qF "distribution" ; then
               # Remove distribution
               rm -f ${file}
            else
               # Check if this file should be compared
               compare=true
               for special in ${not_compared} ; do
                  if printf %s\\n "${suffix}" | grep -qF "${special}" ; then
                     compare=false
                  fi
               done

               if test "${compare}" = "true"; then
                  # Build reference file name, with a special case where the file is mpi-dependent
                  fileref=testref/${test}/test_1-1_${suffix}.nc
                  for special in ${mpi_dependent} ; do
                     if printf %s\\n "${suffix}" | grep -qF "${special}" ; then
                        fileref=testref/${test}/test_${mpi}-1_${suffix}.nc
                     fi
                  done

                  if [ -x "$(command -v nccmp)" ] ; then
                     # Compare files with NCCMP
                     echo -e "Command: nccmp -dfFmqS --threads=${nthreads} -T ${tolerance} ${file} ${fileref}"
                     nccmp -dfFmqS --threads=${nthreads} -T ${tolerance} ${file} ${fileref}
                     exit_code=$?
                     if test "${exit_code}" != "0" ; then
                        echo -e "\e[31mTest failed (nccmp) checking: "${file#testdata/}"\e[0m"
                        status=2
                     fi
                  else
                     echo -e "\e[31mCannot find command: nccmp\e[0m"
                     status=3
                  fi
               fi
            fi
         fi
      done

      # Check number of files
      nfiles=`ls testoutput/${test}/test_${mpi}-${omp}.000000.out 2>/dev/null | wc -l`
      if test "${nfiles}" = "0"; then
         echo -e "\e[31mNo log file to check\e[0m"
         status=4
      fi

      for file in `ls testoutput/${test}/test_${mpi}-${omp}.000000.out 2>/dev/null` ; do
         # Check for stars
         grep -q "*" ${file}
         exit_code=$?
         if test "${exit_code}" = "0" ; then
            echo -e "\e[31mTest failed (stars in output) checking: "${file#testoutput/}"\e[0m"
            status=5
         fi

         # Check for "Close listings" line
         grep -q "\-\-\- Close listings" ${file}
         exit_code=$?
         if test "${exit_code}" != "0" ; then
            echo -e "\e[31mTest failed (no listing closure) checking: "${file#testoutput/}"\e[0m"
            status=6
         fi
      done
   fi

   if test "${compare_type}" = "specific" ; then
      # Specific tests
      testref=$2
      mpiomp=$3

      # Check number of files
      nfiles=`ls testdata/${test}/test_${mpiomp}_*.nc 2>/dev/null | wc -l`
      if test "${nfiles}" = "0"; then
         echo -e "\e[31mNo NetCDF file to check\e[0m"
         status=7
      fi

      for file in `ls testdata/${test}/test_${mpiomp}_*.nc 2>/dev/null` ; do
         if test ! -L ${file}; then
            # Get suffix
            tmp=${file#testdata/${test}/test_${mpiomp}_}
            suffix=${tmp%.nc}

            if printf %s\\n "${suffix}" | grep -qF "distribution" ; then
               # Remove distribution
               rm -f ${file}
            else
               # Check if this file should be compared
               compare=true
               for special in ${not_compared} ; do
                  if printf %s\\n "${suffix}" | grep -qF "${special}" ; then
                     compare=false
                  fi
               done

               if test "${compare}" = "true"; then
                  fileref=testdata/${test}/test_${mpiomp}_${suffix}.nc
                  if [ -x "$(command -v nccmp)" ] ; then
                     # Compare files with NCCMP
                     echo -e "Command: nccmp -dfFmqS --threads=${nthreads} -T ${tolerance} ${file} ${fileref}"
                     nccmp -dfFmqS --threads=${nthreads} -T ${tolerance} ${file} ${fileref}
                     exit_code=$?
                     if test "${exit_code}" != "0" ; then
                        echo -e "\e[31mTest failed (nccmp) checking: "${file#testdata/}"\e[0m"
                        status=8
                     fi
                  else
                     echo -e "\e[31mCannot find command: nccmp\e[0m"
                     status=9
                  fi
               fi
            fi
         fi
      done
   fi
fi

# Exit
exit ${status}
