#!/bin/sh
#----------------------------------------------------------------------
# Shell script: saber_tar
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Parameters
testdir=$1
datadir=$2
list=$3
name=$4

# Loop over files
files=''
while IFS= read -r line
do
   files=${files}' '${line}   
done < ${list}

# Archive
cd ${testdir}/${datadir}
tar -cvzf ${name}.tar.gz ${files}

# Compute MD5 checksum for
cd ${testdir}
md5sum ${datadir}/${name}.tar.gz > ${datadir}/${name}.tar.gz.md5
