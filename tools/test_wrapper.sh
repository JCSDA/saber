#!/bin/bash

# run the executable
echo ""
echo "==============================================================================="
echo "Running test executable"
echo "==============================================================================="

# input parameters
mpiexec_cmd=$1
mpiexec_n=$2
mpi=$3
omp=$4
exename=$5
yamlname=$6
run_file=$7
comparename=$8
ref_file=$9
ftol=${10}
idif=${11}

# Run Test
cmd="${mpiexec_cmd} ${mpiexec_n} ${mpi} ${exename} ${yamlname} ${run_file}"
echo ${cmd}
eval ${cmd}

exit_code=$?
if test "${exit_code}" == "0"; then
    echo -e "Test run succeed"
else
    echo -e "Test run failed with error code: ${exit_code} \n"
    exit ${exit_code}
fi

# Run compare

echo ""
echo "==============================================================================="
echo "Running compare script"
echo "==============================================================================="

cmd="${comparename} ${run_file} ${ref_file} ${ftol} ${idif}"
echo ${cmd}
eval ${cmd}

exit_code=$?
if test "${exit_code}" == "0"; then
   echo -e "Test compare succeed"
else
    echo -e "Test compare failed with error code: ${exit_code}"
    exit ${exit_code}
fi

# Test passed!
echo -e "PASSED"
