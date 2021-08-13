#!/usr/bin/env bash
#----------------------------------------------------------------------
# Bash script: saber_tar_ref
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Parameters
if test "$#" = "0" ; then
   datadir=${HOME}/build/gnu_9.3.0/bundle_RelWithDebInfo/saber/test/testdata
   listdir=${HOME}/code/bundle/saber/test/testlist
else
   datadir=$1
   listdir=$2
fi

# References list
ref_list="
saber_ref_1
saber_ref_2
saber_ref_3
saber_ref_cgal
saber_ref_mpi_1
saber_ref_mpi_2
saber_ref_multi
saber_ref_oops"

# Get git branch
branch=`git rev-parse --abbrev-ref HEAD`

for ref in ${ref_list}; do
   # Loop over files
   files=''
   while IFS= read -r line
   do
      files=${files}' '${line}
   done < ${listdir}/${ref}.txt

   # Go to data directory
   if test ${ref} == "saber_ref_oops"; then
      # Special case for OOPS references
      cd ${datadir}/../testoutput
   else
      cd ${datadir}
   fi

   # Archive
   tar -cvzf ${ref}.tar.gz ${files}

   # Compute MD5 checksum
   md5sum ${ref}.tar.gz > ${ref}.tar.gz.md5

   # Send to S3
   aws s3 cp ${ref}.tar.gz s3://jedi-test-files/saber/${branch}/ --grants read=uri=http://acs.amazonaws.com/groups/global/AllUsers
   aws s3 cp ${ref}.tar.gz.md5 s3://jedi-test-files/saber/${branch}/ --grants read=uri=http://acs.amazonaws.com/groups/global/AllUsers

   # Clean archive and MD5 checksum
   rm -f ${ref}.tar.gz ${ref}.tar.gz.md5
done
