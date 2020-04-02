#!/bin/bash
#----------------------------------------------------------------------
# Bash script: saber_test_wrapper
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Parameters
mpiexec_cmd=$1
mpiexec_n=$2
mpi=$3
omp=$4
exename=$5
yamlname=$6
outputdir=$7
comparename=$8
testname=$9
ctol=${10}
cpu_list=${11}

# Test run
export OMP_NUM_THREADS=${omp}
if test "${cpu_list}" = ""; then
   cmd="${mpiexec_cmd} ${mpiexec_n} ${mpi} ${exename} ${yamlname} ${outputdir}"
else
   cmd="${mpiexec_cmd} ${mpiexec_n} ${mpi} -cpu-list ${cpu_list} ${exename} ${yamlname} ${outputdir}"
fi
eval ${cmd}
exit_code=$?
if test "${exit_code}" == "0" ; then
    echo -e "Test run succeed"
else
    echo -e "Test run failed with error code: ${exit_code} \n"
    exit ${exit_code}
fi

# Test compare
mpixomp=$((mpi*omp))
cmd="${comparename} ${test} ${mpi} ${omp}"
eval ${cmd}
exit_code=$?
if test "${exit_code}" == "0" ; then
   echo -e "Test compare succeed"
else
    echo -e "Test compare failed with error code: ${exit_code}"
    exit ${exit_code}
fi

# Test passed!
echo -e "PASSED"
