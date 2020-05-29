#!/bin/bash
#----------------------------------------------------------------------
# Bash script: saber_parallel
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Check bash version
if test ${BASH_VERSINFO:-0} -lt 4; then
   echo "Bash version lower than 4 does not support associative arrays"
   exit 1
fi

# Check location
PWD=`pwd`
if test "${PWD##*/}" != "test"; then
   echo "This script should be run from \${build_directory}/saber/test"
   exit 2
fi
tmp=${PWD%/test}
if test "${tmp##*/}" != "saber"; then
   echo "This script should be run from \${build_directory}/saber/test"
   exit 2
fi
echo "Working directory: ${PWD}"

# Check number of CPUs
nprocmax=`python -c 'import multiprocessing as mp; print(mp.cpu_count())'`
nproc=$1
if test -z ${nproc}; then
   nproc=${nprocmax}
fi
if test ${nproc} -lt 1; then
   echo "At least one processor is required"
   exit 3
fi
if test ${nproc} -gt ${nprocmax}; then
   echo "Only ${nprocmax} processor(s) are available"
   exit 4
fi
echo "Tests run on ${nproc} logical core(s) (${nprocmax} available)"

# Make temporary directory
mkdir -p saber_ctest_log
rm -fr saber_ctest_log/*

# Initial time
initial_time=`date`
initial_time_sec=`date +%s`
echo "Tests start at ${initial_time}" > saber_ctest_log/execution.log

# List tests
list_get=`ctest -N | grep get_saber | awk '{print $(NF)}'`
list_run=`ctest -N | grep test_bump | grep _run | awk '{print $(NF)}'`
list_compare=`ctest -N | grep test_bump | grep _compare | awk '{print $(NF)}'`
list_qg=`ctest -N | grep test_qg | awk '{print $(NF)}'`

# Tests variables
list_get_array=(${list_get})
ntest_get=${#list_get_array[@]}
stest_get=0
ftest_get=0
list_run_array=(${list_run})
ntest_run=${#list_run_array[@]}
stest_run=0
ftest_run=0
list_compare_array=(${list_compare})
ntest_compare=${#list_compare_array[@]}
stest_compare=0
ftest_compare=0
list_qg_array=(${list_qg})
ntest_qg=${#list_qg_array[@]}
stest_qg=0
ftest_qg=0
ntest_tot=$((ntest_run+ntest_compare+ntest_qg))
itest=0
if test ${ntest_tot} = 0; then
   echo "No test detected, this script should be run from \${build_directory}/saber/test"
   exit 5
fi

# Declare pids array
declare -A pids=()

# Download data
for get in ${list_get}; do
   echo "Handling process ${get}" >> saber_ctest_log/execution.log
   ctest -R ${get} > saber_ctest_log/${get}.log 2> saber_ctest_log/${get}.err &
   pids[${get}]=$!
done

# Wait for all processes to finish
while [ ${#pids[@]} -gt 0 ]; do
   # Sleep
   sleep 0.01

   # Loop over running processes
   for cget in "${!pids[@]}"; do
      # Check if the process is not running
      if ! kill -0 ${pids[${cget}]} 2>/dev/null; then
         # Unset process from the list
         unset pids[${cget}]

         # Check if this process passed
         tmp1=`grep -i "tests passed," saber_ctest_log/${cget}.log`
         tmp2=${tmp1##*tests passed,}
         nerr=${tmp2%% tests failed*}
         if test ${nerr} -eq 0; then
            echo "${cget} passed" >> saber_ctest_log/execution.log
            stest_get=$((stest_get+1))
            echo -e "Download : \033[32mpassed\033[0m ~> ${cget}"
         else
            echo "${cget} failed" >> saber_ctest_log/execution.log
            ftest_get=$((ftest_get+1))
            echo -e "Download : \033[31mfailed\033[0m ~> ${cget}"
         fi
      fi
   done
done

# Declare and initialize cpus array
declare -A cpus=()
for i in $(seq 1 ${nproc}); do
   cpus[$i]=""
done
navail=${nproc}

# Run tests
for run in ${list_run}; do
   echo "Handling process ${run}" >> saber_ctest_log/execution.log

   # Get command and arguments
   ctest --timeout 0.00001 -VV -R ${run} > saber_ctest_log/${run}_fail.log 2>/dev/null
   tmp=`grep "Test command:" saber_ctest_log/${run}_fail.log`
   cmd=`echo ${tmp} | awk '{print $4}'`
   exe=`echo ${tmp} | awk '{print $7,$8,$9}'`
   tmp1=${run%%_run}
   omp=${tmp1##*-}
   tmp2=${tmp1%-*}
   mpi=${tmp2##*_}
   mpixomp=$((mpi*omp))
   if [ ${mpixomp} -gt ${nproc} ]; then
      echo "Too many CPUs are required for process ${run} (${mpixomp} > ${nproc})"
      exit 6
   fi
   rm -f saber_ctest_log/${run}_fail.log

   # Print info
   echo "   At this time, ${navail} CPUs are available and ${mpixomp} are required" >> saber_ctest_log/execution.log

   # Wait for available CPUs
   while [ ${navail} -lt ${mpixomp} ]; do
      # Sleep
      sleep 0.01

      # Loop over running processes
      for crun in "${!pids[@]}"; do
         # Check if the process is not running
         if ! kill -0 ${pids[${crun}]} 2>/dev/null; then
            # Unset process from the list
            unset pids[${crun}]

            # CPUs for this process are now available
            navail=0
            for i in "${!cpus[@]}"; do
               if test "${cpus[${i}]}" = ""; then
                  navail=$((navail+1))
               fi
               if test "${cpus[${i}]}" = "${crun}"; then
                  # This CPU is now available
                  cpus[${i}]=""
                  navail=$((navail+1))
               fi
            done
            echo "   Process ${crun} is done, ${navail} are now are available" >> saber_ctest_log/execution.log

            # Check if this process passed
            err=`wc -l saber_ctest_log/${crun}.err | awk '{print $1}'`
            itest=$((itest+1))
            itest_tot=`printf "%03d" ${itest}`
            if test "${err}" = "0"; then
               # Run passed
               echo "${crun} passed" >> saber_ctest_log/execution.log
               stest_run=$((stest_run+1))
               echo -e "${itest_tot} / ${ntest_tot}: \033[32mpassed\033[0m ~> ${crun}"
            else
               # Run failed
               echo "${crun} failed" >> saber_ctest_log/execution.log
               ftest_run=$((ftest_run+1))
               echo -e "${itest_tot} / ${ntest_tot}: \033[31mfailed\033[0m ~> ${crun}"
            fi
         fi
      done
   done

   # Find available cpus
   cpu_list=""
   ncpus=0
   for i in $(seq 1 ${nproc}); do
      if test "${cpus[${i}]}" = ""; then
         # This CPU is available
         if test ${ncpus} -lt ${mpixomp}; then
            # Use this CPU
            im1=$((i-1))
            cpu_list=${cpu_list}${im1}
            ncpus=$((ncpus+1))
            navail=$((navail-1))

            # Add comma
            if test ${ncpus} -lt ${mpixomp}; then
               cpu_list=${cpu_list}","
            fi

            # Set a task for this CPU
            cpus[${i}]="${run}"
         fi
      fi
   done

   # Fire!
   echo "   Process ${run} is launched on CPUs ${cpu_list}" >> saber_ctest_log/execution.log
   export OMP_NUM_THREADS=${omp}
   full_cmd="${cmd} \"-n\" \"${mpi}\" \"-cpu-list\" \"${cpu_list}\" ${exe} > saber_ctest_log/${run}.log 2> saber_ctest_log/${run}.err &"
   eval ${full_cmd}
   pids[${run}]=$!
done

# Wait for all processes to finish
while [ ${#pids[@]} -gt 0 ]; do
   # Sleep
   sleep 0.01

   # Loop over running processes
   for crun in "${!pids[@]}"; do
      # Check if the process is not running
      if ! kill -0 ${pids[${crun}]} 2>/dev/null; then
         # Unset process from the list
         unset pids[${crun}]

         # Check if this process passed
         err=`wc -l saber_ctest_log/${crun}.err | awk '{print $1}'`
         itest=$((itest+1))
         itest_tot=`printf "%03d" ${itest}`
         if test "${err}" = "0"; then
            # Run passed
            echo "${crun} passed" >> saber_ctest_log/execution.log
            stest_run=$((stest_run+1))
            echo -e "${itest_tot} / ${ntest_tot}: \033[32mpassed\033[0m ~> ${crun}"
         else
            # Run failed
            echo "${crun} failed" >> saber_ctest_log/execution.log
            ftest_run=$((ftest_run+1))
            echo -e "${itest_tot} / ${ntest_tot}: \033[31mfailed\033[0m ~> ${crun}"
         fi
      fi
   done
done

# Re-initialize cpus array
for i in $(seq 1 ${nproc}); do
   cpus[$i]=""
done
navail=${nproc}

# Compare
for compare in ${list_compare}; do
   echo "Handling process ${compare}" >> saber_ctest_log/execution.log

   # Get command and arguments
   ctest --timeout 0.00001 -VV -R ${compare} > saber_ctest_log/${compare}_fail.log 2>/dev/null
   tmp=`grep "Test command:" saber_ctest_log/${compare}_fail.log`
   exe=`echo ${tmp} | awk '{print $4,$5,$6,$7}'`
   mpi=1
   omp=1
   mpixomp=$((mpi*omp))
   rm -f saber_ctest_log/${compare}_fail.log

   # Print info
   echo "   At this time, ${navail} CPUs are available and ${mpixomp} are required" >> saber_ctest_log/execution.log

   # Wait for available CPUs
   while [ ${navail} -lt ${mpixomp} ]; do
      # Sleep
      sleep 0.01

      # Loop over running processes
      for ccompare in "${!pids[@]}"; do
         # Check if the process is not running
         if ! kill -0 ${pids[${ccompare}]} 2>/dev/null; then
            # Unset process from the list
            unset pids[${ccompare}]

            # CPUs for this process are now available
            navail=0
            for i in "${!cpus[@]}"; do
               if test "${cpus[${i}]}" = ""; then
                  navail=$((navail+1))
               fi
               if test "${cpus[${i}]}" = "${ccompare}"; then
                  # This CPU is now available
                  cpus[${i}]=""
                  navail=$((navail+1))
               fi
            done
            echo "   Process ${ccompare} is done, ${navail} are now are available" >> saber_ctest_log/execution.log

            # Check if this process passed
            err=`wc -l saber_ctest_log/${ccompare}.err | awk '{print $1}'`
            itest=$((itest+1))
            itest_tot=`printf "%03d" ${itest}`
            if test "${err}" = "0"; then
               # Compare passed
               echo "${ccompare} passed" >> saber_ctest_log/execution.log
               stest_compare=$((stest_compare+1))
               echo -e "${itest_tot} / ${ntest_tot}: \033[32mpassed\033[0m ~> ${ccompare}"
            else
               # Compare failed
               echo "${ccompare} failed" >> saber_ctest_log/execution.log
               ftest_compare=$((ftest_compare+1))
               echo -e "${itest_tot} / ${ntest_tot}: \033[31mfailed\033[0m ~> ${ccompare}"
            fi
         fi
      done
   done

   # Find available cpus
   cpu_list=""
   ncpus=0
   for i in $(seq 1 ${nproc}); do
      if [ "${cpus[${i}]}" == "" ]; then
         # This CPU is available
         if [ ${ncpus} -lt ${mpixomp} ]; then
            # Use this CPU
            im1=$((i-1))
            cpu_list=${cpu_list}${im1}
            ncpus=$((ncpus+1))
            navail=$((navail-1))

            # Add comma
            if [ ${ncpus} -lt ${mpixomp} ]; then
               cpu_list=${cpu_list}","
            fi

            # Set a task for this CPU
            cpus[${i}]="${compare}"
         fi
      fi
   done

   # Fire!
   echo "   Process ${compare} is launched on CPUs ${cpu_list}" >> saber_ctest_log/execution.log
   export OMP_NUM_THREADS=${omp}
   full_cmd="mpirun \"-n\" \"${mpi}\" \"-cpu-list\" \"${cpu_list}\" ${exe} > saber_ctest_log/${compare}.log 2> saber_ctest_log/${compare}.err &"
   eval ${full_cmd}
   pids[${compare}]=$!
done

# Wait for all processes to finish
while [ ${#pids[@]} -gt 0 ]; do
   # Sleep
   sleep 0.01

   # Loop over running processes
   for ccompare in "${!pids[@]}"; do
      # Check if the process is not running
      if ! kill -0 ${pids[${ccompare}]} 2>/dev/null; then
         # Unset process from the list
         unset pids[${ccompare}]

         # Check if this process passed
         err=`wc -l saber_ctest_log/${ccompare}.err | awk '{print $1}'`
         itest=$((itest+1))
         itest_tot=`printf "%03d" ${itest}`
         if test "${err}" = "0"; then
            # Compare passed
            echo "${ccompare} passed" >> saber_ctest_log/execution.log
            stest_compare=$((stest_compare+1))
            echo -e "${itest_tot} / ${ntest_tot}: \033[32mpassed\033[0m ~> ${ccompare}"
         else
            # Compare failed
            echo "${ccompare} failed" >> saber_ctest_log/execution.log
            ftest_compare=$((ftest_compare+1))
            echo -e "${itest_tot} / ${ntest_tot}: \033[31mfailed\033[0m ~> ${ccompare}"
         fi
      fi
   done
done

# Re-initialize cpus array
for i in $(seq 1 ${nproc}); do
   cpus[$i]=""
done
navail=${nproc}

# QG
for qg in ${list_qg}; do
   echo "Handling process ${qg}" >> saber_ctest_log/execution.log

   # Get command and arguments
   ctest -VV -R ${qg} > saber_ctest_log/${qg}.log 2> saber_ctest_log/${qg}.err

   # Check if this process passed
   err=`grep PASSED saber_ctest_log/${qg}.log | awk '{print $2}'`
   itest=$((itest+1))
   itest_tot=`printf "%03d" ${itest}`
   if test "${err}" = "PASSED"; then
      # QG passed
      echo "${qg} passed" >> saber_ctest_log/execution.log
      stest_qg=$((stest_qg+1))
      echo -e "${itest_tot} / ${ntest_tot}: \033[32mpassed\033[0m ~> ${qg}"
   else
      # QG failed
      echo "${qg} failed" >> saber_ctest_log/execution.log
      ftest_qg=$((ftest_qg+1))
      echo -e "${itest_tot} / ${ntest_tot}: \033[31mfailed\033[0m ~> ${qg}"
   fi
done

# Final time
final_time=`date`
final_time_sec=`date +%s`
echo "Tests finished at ${final_time}" >> saber_ctest_log/execution.log

# Elapsed time
elapsed_time=$((final_time_sec-initial_time_sec))
echo "Elapsed time: ${elapsed_time} sec." >> saber_ctest_log/execution.log
echo -e ""
echo -e "Elapsed time: \033[36m${elapsed_time}\033[0m sec."

# Grep failed tests
failed_tests=`grep failed saber_ctest_log/execution.log | awk '{print $1}'`
ftest=$((ftest_get+ftest_run+ftest_compare+ftest_qg))
if test "${ftest}" -gt "0"; then
   echo -e ""
   echo -e "Failed tests:"
   for failed_test in ${failed_tests}; do
      echo -e "  \033[31m${failed_test}\033[0m"
   done
fi

# Print summary
echo -e ""
echo -e "Summary:"
# Download tests
stest=`printf "%03d" $((stest_get))`
ftest=`printf "%03d" $((ftest_get))`
echo -e "  Download: \033[32m${stest}\033[0m tests passed and \033[31m${ftest}\033[0m failed"

# Run tests
stest=`printf "%03d" $((stest_run))`
ftest=`printf "%03d" $((ftest_run))`
echo -e "  Run:      \033[32m${stest}\033[0m tests passed and \033[31m${ftest}\033[0m failed"

# Compare tests
stest=`printf "%03d" $((stest_compare))`
ftest=`printf "%03d" $((ftest_compare))`
echo -e "  Compare:  \033[32m${stest}\033[0m tests passed and \033[31m${ftest}\033[0m failed"

# QG tests
stest=`printf "%03d" $((stest_qg))`
ftest=`printf "%03d" $((ftest_qg))`
echo -e "  QG:       \033[32m${stest}\033[0m tests passed and \033[31m${ftest}\033[0m failed"

exit 0
