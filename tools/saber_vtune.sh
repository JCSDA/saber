#!/bin/bash
# This bash script can be used to assess code performance with Intel VTune.
# First, compile the code with a recent Intel compiler.
# Second, run "ctest -N" to determine the index of the test you wish to evaluate.
# Finally, update this script with your own build directory ('build_dir' variable) and run it,
# specifyng the type of VTune diagnostic (hotspots, memory-access or hpc-performance) and the
# test index as arguments.
# The optional argument 'final_tar' can be used to produce a tar file with the final results instead of
# opening them with vtune-gui.
#
# Usage: saber_vtune.sh diagnostic test_index [final_tar]

if [[ $# -ne 2 ]]; then
   echo "Error: wrong number of arguments"
   echo "Usage: saber_vtune.sh diagnostic test_index [final_tar]"
   echo "- diagnostic can be: hotspots, memory-access or hpc-performance"
   echo "- type \"ctest -N\" to list the test indices"
   echo "- final_tar is optional and can be true or false (default is false), if set to true, the final results are gathered into a tar file"
   exit 1
fi

# Input paramters
diagnostic=$1
test_index=$2
if [[ $# -eq 1 ]]; then
   final_tar=$3
else
   final_tar=false
fi

# Hard-coded parameters
build_dir=${HOME}/build/intel_2021.5.0/bundle_RelWithDebInfo

# Go to SABER build directory
cd ${build_dir}/saber/test

indices=`ctest -N | grep saber_test_ | gawk '{print $2}'`
for index in ${indices}; do
   name=`ctest -N | grep "${index}" | gawk '{print $3}'`
   index=${index##\#}
   index=${index%%:*}

   if [[ ${index} -eq ${test_index} ]]; then
      valid=true
   else
      valid=false
   fi

   if [[ ${valid} == "true" ]]; then
      echo "Test #"${index}": "${name}
      rm -fr vtune.${name}.${diagnostic}.${HOSTNAME}
      ctest --timeout 0.00001 -VV -I ${index},${index} > tmp.log.out 2>/dev/null
      tmp=`grep "Test command:" tmp.log.out`
      exe=`echo ${tmp} | awk '{print $7,$8,$9}'`
      omp=${name##*-}
      tmp2=${name%-*}
      mpi=${tmp2##*_}
      full_cmd="mpirun \"-n\" \"${mpi}\" -l vtune -quiet -collect ${diagnostic} -trace-mpi -result-dir vtune.${name}.${diagnostic} ${exe} > tmp.log.out 2>/dev/null"
      echo ${full_cmd}
      eval ${full_cmd}
      rm -f tmp.log.out
      if [[ ${final_tar} == "true" ]]; then
         # Tar results
         tar -cf vtune.${name}.${diagnostic}.${HOSTNAME}.tar vtune.${name}.${diagnostic}.${HOSTNAME}
      elif [[ ${final_tar} == "false" ]]; then
         # Visualize with vtune-gui
         vtune-gui vtune.${name}.${diagnostic}.${HOSTNAME} &
      else
         # Wrong final_tar argument
         echo "Error: wrong final_tar argument, should be true or false"
         exit 2
      fi
   fi
done

exit 0
