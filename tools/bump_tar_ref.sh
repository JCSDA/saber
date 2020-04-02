#!/bin/bash
#----------------------------------------------------------------------
# Bash script: bump_tar_ref
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Parameters
datadir=$1
listdir=$2

# References list
ref_list="
bump_ref_1
bump_ref_2
bump_ref_3
bump_ref_mpi_1
bump_ref_mpi_2
bump_ref_mpi_3
bump_ref_quad"

# Get git branch
branch=`git rev-parse --abbrev-ref HEAD`

# Get initial pwd
ipwd=`pwd`

for ref in ${ref_list}; do
   # Loop over files
   files=''
   while IFS= read -r line
   do
      files=${files}' '${line}   
   done < ${listdir}/${ref}.txt

   # Archive
   cd ${datadir}
   tar -cvzf ${ref}.tar.gz ${files}

   # Compute MD5 checksum
   md5sum ${ref}.tar.gz > ${ref}.tar.gz.md5

   # Send to S3
   aws s3 cp ${ref}.tar.gz s3://jedi-test-files/saber/${branch}/ --grants read=uri=http://acs.amazonaws.com/groups/global/AllUsers
   aws s3 cp ${ref}.tar.gz.md5 s3://jedi-test-files/saber/${branch}/ --grants read=uri=http://acs.amazonaws.com/groups/global/AllUsers

   # Clean archive and MD5 checksum
   rm -f ${ref}.tar.gz ${ref}.tar.gz.md5

   # Back to initial directory
   cd ${ipwd}
done
