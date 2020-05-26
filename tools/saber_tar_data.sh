#!/bin/bash
#----------------------------------------------------------------------
# Bash script: saber_tar_data
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Parameters
datadir=$1
listdir=$2

# Data list
data_list="
saber_data
saber_data_mpi
saber_data_omp
saber_data_oops"

# Get git branch
branch=`git rev-parse --abbrev-ref HEAD`

# Get initial pwd
ipwd=`pwd`

for data in ${data_list}; do
   # Loop over files
   files=''
   while IFS= read -r line
   do
      line_tmp="$(echo ${line} | sed 's/read/write/g')"
      if test "${line}" != "${line_tmp}"; then
         mv -f ${datadir}/${line} ${datadir}/${line}_old
         cp -f ${datadir}/${line_tmp} ${datadir}/${line}
      fi
      files=${files}' '${line}
   done < ${listdir}/${data}.txt

   # Archive
   cd ${datadir}
   tar -hcvzf ${data}.tar.gz ${files}

   # Compute MD5 checksum
   md5sum ${data}.tar.gz > ${data}.tar.gz.md5

   # Send to S3
   aws s3 cp ${data}.tar.gz s3://jedi-test-files/saber/${branch}/ --grants read=uri=http://acs.amazonaws.com/groups/global/AllUsers
   aws s3 cp ${data}.tar.gz.md5 s3://jedi-test-files/saber/${branch}/ --grants read=uri=http://acs.amazonaws.com/groups/global/AllUsers

   # Clean archive and MD5 checksum
   rm -f ${data}.tar.gz ${data}.tar.gz.md5

   # Back to initial directory
   cd ${ipwd}

   # Loop over files
   files=''
   while IFS= read -r line
   do
      line_tmp="$(echo ${line} | sed 's/read/write/g')"
      if test "${line}" != "${line_tmp}"; then
         rm -f ${datadir}/${line}
         mv -f ${datadir}/${line}_old ${datadir}/${line}
      fi
   done < ${listdir}/${data}.txt
done
