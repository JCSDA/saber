#!/usr/bin/env bash
#----------------------------------------------------------------------
# Bash script: saber_tar_ref
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Default parameters
datadir=${HOME}/build/gnu_9.3.0/bundle_RelWithDebInfo/saber/test/testdata
listdir=${HOME}/code/bundle/saber/test/testlist
branch=`git rev-parse --abbrev-ref HEAD`
upload_on_s3=0

# Input parameters
if test "$#" -gt "0" ; then
   datadir=$1
   listdir=$2
   branch=$3
   upload_on_s3=$4
fi

# saber-data/testref directory
saber_data=${listdir}/../../../saber-data/testref

if test ${upload_on_s3} == "1"; then
   # Tier 1 to 3 tests
   for tier in $(seq 1 3); do
      # Loop over files
      files=''
      while IFS= read -r line
      do
         # OOPS-specific tests
         if test ${line#qg_} != ${line}; then
            ln -sf ${datadir}/../testoutput/${line} ${datadir}/${line}
         fi

         # Add file
         files=${files}' '${line}
      done < ${listdir}/saber_ref_tier${tier}.txt

      # Go to data directory
      cd ${datadir}

      # Archive
      tar -hcvzf saber_ref_tier${tier}.tar.gz ${files}

      # Compute MD5 checksum
      md5sum saber_ref_tier${tier}.tar.gz > saber_ref_tier${tier}.tar.gz.md5

      # Send to S3
      aws s3 cp saber_ref_tier${tier}.tar.gz s3://jedi-test-files/saber/${branch}/ --grants read=uri=http://acs.amazonaws.com/groups/global/AllUsers
      aws s3 cp saber_ref_tier${tier}.tar.gz.md5 s3://jedi-test-files/saber/${branch}/ --grants read=uri=http://acs.amazonaws.com/groups/global/AllUsers

      # Clean archive and MD5 checksum
      rm -f saber_ref_tier${tier}.tar.gz saber_ref_tier${tier}.tar.gz.md5
   done
else
   # Tier 1 to 3 tests
   for tier in $(seq 1 3); do
      # Loop over files
      while IFS= read -r line
      do
         # OOPS-specific tests
         if test ${line#qg_} != ${line}; then
            ln -sf ${datadir}/../testoutput/${line} ${datadir}/${line}
         fi

         # Copy files to saber-data/testref
         mkdir -p `dirname ${saber_data}/${line}`
         echo "Copy ${line}"
         cp -f ${datadir}/${line} ${saber_data}/${line}
      done < ${listdir}/saber_ref_tier${tier}.txt
   done
   cd ${saber_data}
   git status
fi
