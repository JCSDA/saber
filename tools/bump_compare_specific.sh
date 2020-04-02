#!/bin/bash
#----------------------------------------------------------------------
# Bash script: bump_compare_specific
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Parameters
test1=$1
test2=$2
suffix=$3

# Initialize exit status
status=0

# NCCMP Parameters
tolerance=1.e-5

# Build file names
file1=testref/${test1}/test_1-1_${suffix}.nc
file2=testref/${test2}/test_1-1_${suffix}.nc

# Compare files with NCCMP
if test -x "$(command -v nccmp)"; then
   echo "Command: nccmp -dfFmqS -T ${tolerance} ${file1} ${file2}"
   nccmp -dfFmqS -T ${tolerance} ${file1} ${file2}
   exit_code=$?
   if test "${exit_code}" != "0"; then
      echo "\e[31mTest failed (nccmp) checking: "${file#testdata/}"\e[0m"
      status=1
      exit ${status}
   fi
else
   echo "\e[31mCannot find command: nccmp\e[0m"
   status=2
   exit ${status}
fi

# Test passed!
echo -e "PASSED"

# Exit
exit ${status}
