#!/bin/sh
#----------------------------------------------------------------------
# Shell script: saber_tar
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Parameters
dir=$1
list=$2
name=$3

# Go to data directory
cd ${dir}

# Loop over files
files=''
while IFS= read -r line
do
   files=${files}' '${line}   
done < ${list}

# Archive
tar -cvzf ${name}.tar.gz ${files}
