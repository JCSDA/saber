#!/usr/bin/env bash
#----------------------------------------------------------------------
# Bash script: saber_setup
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Get arguments
TESTDATA_DIR=$1
TESTREF_DIR=$2
CMAKE_CURRENT_BINARY_DIR=$3
CMAKE_BINARY_DIR=$4

# Clear working directories
rm -fr ${CMAKE_CURRENT_BINARY_DIR}/testdata
rm -fr ${CMAKE_CURRENT_BINARY_DIR}/testinput
rm -fr ${CMAKE_CURRENT_BINARY_DIR}/testoutput
rm -fr ${CMAKE_CURRENT_BINARY_DIR}/testref

# Make working directories
mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/testdata
mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/testinput
mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/testoutput
mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/testref
while IFS= read -r line; do
   mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/testdata/${line}
   mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/testoutput/${line}
   mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/testref/${line}
done < ${CMAKE_BINARY_DIR}/bin/saber_testdir

# Link data files
while IFS= read -r line; do
   ln -sf ${TESTDATA_DIR}/${line} ${CMAKE_CURRENT_BINARY_DIR}/testdata/${line}
done < ${CMAKE_BINARY_DIR}/bin/saber_testdata

# Link reference file
while IFS= read -r line; do
   ln -sf ${TESTREF_DIR}/${line} ${CMAKE_CURRENT_BINARY_DIR}/testref/${line}
done < ${CMAKE_BINARY_DIR}/bin/saber_testref
