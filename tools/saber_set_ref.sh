#!/bin/sh
#----------------------------------------------------------------------
# Shell script: saber_set_ref
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Parameters
testdata=$1
listdir=$2

# Special suffixes list
special_list="mom lct_cor nicas normality obs sampling_grids vbal"

for tier in $(seq 1 3); do
   # Initialize lists
   rm -f ${listdir}/saber_ref_${tier}.txt
   rm -f ${listdir}/saber_ref_mpi_${tier}.txt

   # Get list of tests
   while IFS= read -r saber_test
   do
      # BUMP tests
      if test "${saber_test%%_*}" = "bump"; then
         # Copy 1-1 files
         for file in `ls ${testdata}/${saber_test}/test_1-1_*.nc`; do
            if test ! -L ${file}; then
               echo ${saber_test}/"$(basename -- ${file})" >> ${listdir}/saber_ref_${tier}.txt
            fi
         done

         # Copy 2-1 special files
         for special in ${special_list} ; do
            if ls ${testdata}/${saber_test}/test_2-1_${special}*.nc 1> /dev/null 2>&1; then
               for file in `ls ${testdata}/${saber_test}/test_2-1_${special}*.nc`; do
                  if test ! -L ${file}; then
                     echo ${saber_test}/"$(basename -- $file)" >> ${listdir}/saber_ref_mpi_${tier}.txt
                  fi
               done
            fi
         done
      fi
   done < ${listdir}/saber_test_${tier}.txt
done

# Quad-core tests

# Initialize list
rm -f ${listdir}/saber_ref_quad.txt

# Get list of tests
while IFS= read -r saber_test
do
   # BUMP tests
   if test "${saber_test%%_*}" = "bump"; then
      # Copy 4-1 special files
      for special in ${special_list} ; do
         if ls ${testdata}/${saber_test}/test_4-1_${special}*.nc 1> /dev/null 2>&1; then
            for file in `ls ${testdata}/${saber_test}/test_4-1_${special}*.nc`; do
               if test ! -L ${file}; then
                  echo ${saber_test}/"$(basename -- $file)" >> ${listdir}/saber_ref_quad.txt
               fi
            done
         fi
      done
   fi
done < ${listdir}/saber_test_quad.txt

# OOPS tests

# Initialize list
rm -f ${listdir}/saber_ref_oops.txt

# Get list of tests
while IFS= read -r oops_test
do
   # OOPS tests
   echo ${oops_test}/"test.log.out" >> ${listdir}/saber_ref_oops.txt
done < ${listdir}/saber_test_oops.txt
