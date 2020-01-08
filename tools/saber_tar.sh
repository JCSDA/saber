#!/bin/sh
#----------------------------------------------------------------------
# Shell script: saber_tar
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Parameters
datadir=$1
listdir=$2
name=$3

# Get git branch
branch=`git rev-parse --abbrev-ref HEAD`

# Loop over files
files=''
while IFS= read -r line
do
   files=${files}' '${line}   
done < ${listdir}/${name}.txt

# Archive
cd ${datadir}
tar -cvzf ${name}.tar.gz ${files}

# Compute MD5 checksum
md5sum ${name}.tar.gz > ${name}.tar.gz.md5

# Send to S3
aws s3 cp ${name}.tar.gz s3://jedi-test-files/saber/${branch}/ --grants read=uri=http://acs.amazonaws.com/groups/global/AllUsers
aws s3 cp ${name}.tar.gz.md5 s3://jedi-test-files/saber/${branch}/ --grants read=uri=http://acs.amazonaws.com/groups/global/AllUsers
