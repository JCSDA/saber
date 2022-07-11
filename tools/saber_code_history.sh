#!/bin/bash
#----------------------------------------------------------------------
# Bash shell script: saber_code_history
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Root directory
rootdir=`pwd`/..

# Temporary directory
tmpdir=${HOME}/tmp

# Extract old commits data
cd ${rootdir}/docs/history
tar -xzf history.tar.gz

# Clone SABER
cd ${tmpdir}
rm -fr saber
git clone https://github.com/jcsda-internal/saber
cd saber

# Get the list of merge commits
merge_list=()
for commit in `git log --pretty=oneline | gawk '{print $1}'`; do
   merge_commit=`git find-merge ${commit}`
   if [[ ! " ${merge_list[*]} " =~ " ${merge_commit} " ]]; then
      merge_list+=("${merge_commit}")
   fi
done

# Loop over all commits
for commit in `git log --pretty=oneline | gawk '{print $1}'`; do
   # Get commit date
   commit_date=`git show -s --format=%cd --date=format:"%Y%m%d%H%M%S" ${commit}`

   # Check if it is a merge commit
   if [[ "${merge_list[*]}" =~ "${commit}" ]]; then
      echo ${commit_date} ${commit}

      # Check if the commit is already available in the archive
      if test ! -f "${rootdir}/docs/history/${commit_date}.csv"; then
         git checkout -q ${commit}
         cloc --force-lang="Fortran 90",fypp --unix --skip-uniqueness --quiet --csv --out=${rootdir}/docs/history/${commit_date}.csv .
      fi

      # Get number of lines for Fortran, C++ and Python code
      nfortran=`grep -si Fortran ${rootdir}/docs/history/${commit_date}.csv | gawk -F',' '{print $5}' | gawk '{s+=$1} END {print s}'`
      nfortran=${nfortran:-0}
      ncpp=`grep -si C++ ${rootdir}/docs/history/${commit_date}.csv | gawk -F',' '{print $5}' | gawk '{s+=$1} END {print s}'`
      ncpp=${ncpp:-0}
      npython=`grep -si Python ${rootdir}/docs/history/${commit_date}.csv | gawk -F',' '{print $5}' | gawk '{s+=$1} END {print s}'`
      npython=${npython:-0}

      # Fill database
      commit_date_detail=`git show -s --format=%cd --date=format:"%Y %m %d %H %M %S" ${commit}`
      echo ${commit_date_detail} ${nfortran} ${ncpp} ${npython} >> ${rootdir}/docs/history/history.txt
   fi
done

# Archive old commits data
cd ${rootdir}/docs/history
rm -f history.tar.gz
tar -czf history.tar.gz *.csv
rm -f *.csv

# Run python script to plot history
python3 ${rootdir}/tools/saber_code_history.py

# Resize image
mogrify -resize 25% ${rootdir}/docs/history/history.jpg

# Cleaning
rm -f history.txt
rm -fr ${tmpdir}/saber
