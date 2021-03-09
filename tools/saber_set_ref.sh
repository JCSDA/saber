#!/usr/bin/env bash
#----------------------------------------------------------------------
# Bash script: saber_set_ref
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Parameters
if test "$#" = "0" ; then
   datadir=${HOME}/build/gnu_9.3.0/bundle/saber/test/testdata
   listdir=${HOME}/code/bundle/saber/test/testlist
else
   datadir=$1
   listdir=$2
fi

# Special suffixes list for BUMP
special_list="mom lct_cor nicas normality sampling_grids obs vbal"

# Multi-processor tests
multi_list=$(seq 4 2 12)

# Tier 1 to 3 tests for BUMP
for tier in $(seq 1 3); do
   # Remove lists
   rm -f ${listdir}/saber_ref_${tier}.txt
   rm -f ${listdir}/saber_ref_mpi_${tier}.txt

   # Loop over tests
   while IFS= read -r bump_test
   do
      # Copy 1-1 files
      for file in `ls ${datadir}/${bump_test}/test_1-1_*.nc`; do
         if test ! -L ${file}; then
            echo ${bump_test}/"$(basename -- ${file})" >> ${listdir}/saber_ref_${tier}.txt
         fi
      done

      # Copy 2-1 special files
      for special in ${special_list}; do
         if ls ${datadir}/${bump_test}/test_2-1_${special}*.nc 1> /dev/null 2>&1; then
            for file in `ls ${datadir}/${bump_test}/test_2-1_${special}*.nc`; do
               if test ! -L ${file}; then
                  echo ${bump_test}/"$(basename -- $file)" >> ${listdir}/saber_ref_mpi_${tier}.txt
               fi
            done
         fi
      done
   done < ${listdir}/saber_test_${tier}.txt
done

# CGAL-specific tests for BUMP

# Remove lists
rm -f ${listdir}/saber_ref_cgal.txt
rm -f ${listdir}/saber_ref_mpi_cgal.txt

# Loop over tests
while IFS= read -r bump_test
do
   # Copy 1-1 files
   for file in `ls ${datadir}/${bump_test}/test_1-1_*.nc`; do
      if test ! -L ${file}; then
         echo ${bump_test}/"$(basename -- ${file})" >> ${listdir}/saber_ref_cgal.txt
      fi
   done

   # Copy 2-1 special files
   for special in ${special_list}; do
      if ls ${datadir}/${bump_test}/test_2-1_${special}*.nc 1> /dev/null 2>&1; then
         for file in `ls ${datadir}/${bump_test}/test_2-1_${special}*.nc`; do
            if test ! -L ${file}; then
               echo ${bump_test}/"$(basename -- $file)" >> ${listdir}/saber_ref_mpi_cgal.txt
            fi
         done
      fi
   done
done < ${listdir}/saber_test_cgal.txt

# Multi-core tests for BUMP

# Remove list
rm -f ${listdir}/saber_ref_multi.txt

# Loop over tests
while IFS= read -r bump_test
do
   for multi in ${multi_list}; do
      # Copy N-1 special files
      for special in ${special_list}; do
         if ls ${datadir}/${bump_test}/test_${multi}-1_${special}*.nc 1> /dev/null 2>&1; then
            for file in `ls ${datadir}/${bump_test}/test_${multi}-1_${special}*.nc`; do
               if test ! -L ${file}; then
                  echo ${bump_test}/"$(basename -- $file)" >> ${listdir}/saber_ref_multi.txt
               fi
            done
         fi
      done
   done
done < ${listdir}/saber_test_multi.txt

# Tests for OOPS

# Remove list
rm -f ${listdir}/saber_ref_oops.txt

# Loop over tests
while IFS= read -r oops_test
do
   echo ${oops_test}/"test.log.out" >> ${listdir}/saber_ref_oops.txt
done < ${listdir}/saber_test_oops.txt
