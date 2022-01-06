#!/usr/bin/env bash
#----------------------------------------------------------------------
# Bash script: saber_set_ref
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Parameters
if test "$#" = "0" ; then
   datadir=${HOME}/build/gnu_10.3.0/bundle_RelWithDebInfo/saber/test/testdata
   listdir=${HOME}/code/bundle/saber/test/testlist
else
   datadir=$1
   listdir=$2
fi

# Special suffixes list
mpi_dependent="mom lct_cor nicas_grids nicas_local normality sampling_local sampling_grids vbal_cov_local vbal_local"

# Multi-processor tests
multi_list=$(seq 4 2 12)

# Tier 1 to 3 tests
for tier in $(seq 1 3); do
   # Remove lists
   rm -f ${listdir}/saber_ref_tier${tier}.txt

   # Loop over tests
   while IFS= read -r bump_test
   do
      # Copy 1-1 files
      for file in `ls ${datadir}/${bump_test}/test_1-1_*.nc`; do
         if test ! -L ${file}; then
            echo ${bump_test}/"$(basename -- ${file})" >> ${listdir}/saber_ref_tier${tier}.txt
         fi
      done

      # Copy 2-1 special files
      for special in ${mpi_dependent}; do
         if ls ${datadir}/${bump_test}/test_2-1_${special}*.nc 1> /dev/null 2>&1; then
            for file in `ls ${datadir}/${bump_test}/test_2-1_${special}*.nc`; do
               if test ! -L ${file}; then
                  echo ${bump_test}/"$(basename -- $file)" >> ${listdir}/saber_ref_tier${tier}.txt
               fi
            done
         fi
      done
   done < ${listdir}/saber_test_tier${tier}.txt
done

# Tier 1 CGAL-specific tests

# Loop over tests
while IFS= read -r bump_test
do
   # Copy 1-1 files
   for file in `ls ${datadir}/${bump_test}/test_1-1_*.nc`; do
      if test ! -L ${file}; then
         echo ${bump_test}/"$(basename -- ${file})" >> ${listdir}/saber_ref_tier1.txt
      fi
   done
done < ${listdir}/saber_test_tier1-cgal.txt

# Tier 1 multicore tests

# Loop over tests
while IFS= read -r bump_test
do
   for multi in ${multi_list}; do
      # Copy N-1 special files
      for special in ${mpi_dependent}; do
         if ls ${datadir}/${bump_test}/test_${multi}-1_${special}*.nc 1> /dev/null 2>&1; then
            for file in `ls ${datadir}/${bump_test}/test_${multi}-1_${special}*.nc`; do
               if test ! -L ${file}; then
                  echo ${bump_test}/"$(basename -- $file)" >> ${listdir}/saber_ref_tier1.txt
               fi
            done
         fi
      done
   done
done < ${listdir}/saber_test_tier1-multi.txt

# Tier 1 OOPS-specific tests

# Loop over tests
while IFS= read -r oops_test
do
   echo ${oops_test}/"test.log.out" >> ${listdir}/saber_ref_tier1.txt
done < ${listdir}/saber_test_tier1-oops.txt
