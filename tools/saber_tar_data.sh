#!/usr/bin/env bash
#----------------------------------------------------------------------
# Bash script: saber_tar_data
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Default paramters
datadir=${HOME}/build/gnu_10.3.0/bundle_RelWithDebInfo/saber/test/testdata
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

# saber-data/testdata directory
saber_data=${listdir}/../../../saber-data/testdata

# Get initial pwd
ipwd=`pwd`

if test ${upload_on_s3} == "1"; then
   # Loop over files
   files=''
   while IFS= read -r line
   do
      # Read specific case
      line_tmp="$(echo ${line} | sed 's/read/write/g')"
      if test "${line}" != "${line_tmp}"; then
         mv -f ${datadir}/${line} ${datadir}/${line}_old
         cp -f ${datadir}/${line_tmp} ${datadir}/${line}
      fi

      # Add file
      files=${files}' '${line}
   done < ${listdir}/saber_data.txt

   # Archive
   cd ${datadir}
   tar -hcvzf saber_data.tar.gz ${files}

   # Compute MD5 checksum
   md5sum saber_data.tar.gz > saber_data.tar.gz.md5

   # Send to S3
   aws s3 cp saber_data.tar.gz s3://jedi-test-files/saber/${branch}/ --grants read=uri=http://acs.amazonaws.com/groups/global/AllUsers
   aws s3 cp saber_data.tar.gz.md5 s3://jedi-test-files/saber/${branch}/ --grants read=uri=http://acs.amazonaws.com/groups/global/AllUsers

   # Clean archive and MD5 checksum
   rm -f saber_data.tar.gz saber_data.tar.gz.md5
else
   # Loop over files
   while IFS= read -r line
   do
      # Read specific case
      line_tmp="$(echo ${line} | sed 's/read/write/g')"
      if test "${line}" != "${line_tmp}"; then
         mv -f ${datadir}/${line} ${datadir}/${line}_old
         cp -f ${datadir}/${line_tmp} ${datadir}/${line}
      fi

      # Copy files to saber-data/testdata if files are different
      if ! [ "${datadir}/${line}" -ef "${saber_data}/${line}" ]; then
         mkdir -p `dirname ${saber_data}/${line}`
         echo "Copy ${line}"
         cp -f ${datadir}/${line} ${saber_data}/${line}
      fi
   done < ${listdir}/saber_data.txt
   cd ${saber_data}
   git status
fi

# Back to initial directory
cd ${ipwd}

# Reset files
while IFS= read -r line
do
   # Read specific case
   line_tmp="$(echo ${line} | sed 's/read/write/g')"
   if test "${line}" != "${line_tmp}"; then
      rm -f ${datadir}/${line}
      mv -f ${datadir}/${line}_old ${datadir}/${line}
   fi
done < ${listdir}/saber_data.txt
