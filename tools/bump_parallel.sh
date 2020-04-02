#!/bin/bash
#----------------------------------------------------------------------
# Bash script: bump_parallel
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Check number of CPUs
nproc=$1
if test -z ${nproc}; then
   echo "Specify a number of processors"
   exit 1
fi
if test ${nproc} -lt 1; then
   echo "At least one processor is required"
   exit 2
fi
nprocmax=`nproc --all`
if test ${nproc} -gt ${nprocmax}; then
   echo "Max number of processors is "${nprocmax}
   exit 3
fi

# Make temporary directory
mkdir -p bump_ctest_log
rm -fr bump_ctest_log/*

# Initial time
initial_time=`date --rfc-3339=seconds`
initial_time_sec=`date +%s`
echo "Tests start at ${initial_time}" > bump_ctest_log/execution.log

# Lists
list_get=`ctest -N | grep get_saber | awk '{print $(NF)}'`
list_test=`ctest -N | grep test_bump | awk '{print $(NF)}'`
list_compare=`ctest -N | grep compare_bump | awk '{print $(NF)}'`
list_test_array=(${list_test})
ntest=${#list_test_array[@]}

# Tests variables
stest=0
ftest=0

# Declare pids array
declare -A pids=()

# Download data
for get in ${list_get}; do
   echo "Handling process ${get}" >> bump_ctest_log/execution.log
   ctest -R ${get} > bump_ctest_log/${get}.log &
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

         # Check if this process succeed
         tmp1=`grep -i "tests passed," bump_ctest_log/${cget}.log`
         tmp2=${tmp1##*tests passed,}
         nerr=${tmp2%% tests failed*}
         if test ${nerr} -eq 0; then
            echo -e "Download : \033[32mPassed\033[0m ~> ${cget}"
         else
            echo -e "Download : \033[31mFailed\033[0m ~> ${cget}"
         fi
      fi
   done
done
echo "Data download successful" >> bump_ctest_log/execution.log

# Declare and initialize cpus array
declare -A cpus=()
for i in $(seq 1 ${nproc}); do
   cpus[$i]=""
done
navail=${nproc}

# Run tests
for run in ${list_test}; do
   echo "Handling process ${run}" >> bump_ctest_log/execution.log

   # Get command and arguments
   ctest --timeout 0.00001 -VV -R ${run} > bump_ctest_log/${run}_fail.log 2>/dev/null
   tmp=`grep "Test command:" bump_ctest_log/${run}_fail.log`
   cmd=${tmp##*Test command:}
   mpi=`echo ${cmd} | gawk '{print $4}' | tr -d '"'`
   omp=`echo ${cmd} | gawk '{print $5}' | tr -d '"'`
   mpixomp=$((mpi*omp))
   if test ${mpixomp} -gt ${nproc}; then
      echo "Too many CPUs are required for process ${run} (${mpixomp} > ${nproc})"
      exit 1
   fi

   # Print info
   echo "   At this time, ${navail} CPUs are available and ${mpixomp} are required" >> bump_ctest_log/execution.log

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
            echo "   Process ${crun} is done, ${navail} are now are available" >> bump_ctest_log/execution.log

            # Check if this process succeed
            msg=`tail -n 1 bump_ctest_log/${crun}.log | gawk '{print $1}'`
            if test "${msg}" = "PASSED"; then
               # Run succeed
               stest=$((stest+1))
               itest=`printf "%03d" $((stest+ftest))`
               echo -e "${itest} / ${ntest}: \033[32mPassed\033[0m ~> ${crun}"
            else
               # Run failed
               ftest=$((ftest+1))
               itest=`printf "%03d" $((stest+ftest))`
               echo -e "${itest} / ${ntest}: \033[31mFailed\033[0m ~> ${crun}"
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
   echo "   Process ${run} is launched on CPUs ${cpu_list}" >> bump_ctest_log/execution.log
   full_cmd="${cmd} ${cpu_list} > bump_ctest_log/${run}.log 2> bump_ctest_log/${run}.err &"
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

         # Check if this process succeed
         msg=`tail -n 1 bump_ctest_log/${crun}.log | gawk '{print $1}'`
         if test "${msg}" = "PASSED"; then
            # Run succeed
            stest=$((stest+1))
            itest=`printf "%03d" $((stest+ftest))`
            echo -e "${itest} / ${ntest}: \033[32mPassed\033[0m ~> ${crun}"
         else
            # Run failed
            ftest=$((ftest+1))
            itest=`printf "%03d" $((stest+ftest))`
            echo -e "${itest} / ${ntest}: \033[31mFailed\033[0m ~> ${crun}"
         fi
      fi
   done
done
echo "All runs done" >> bump_ctest_log/execution.log

# Compare tests
for run in ${list_compare}; do
   echo "Handling process ${run}" >> bump_ctest_log/execution.log

   # Run test
   ctest -VV -R ${run} > bump_ctest_log/${run}.log 2>/dev/null

   # Check if this process succeed
   msg=`tail -n 1 bump_ctest_log/${crun}.log | gawk '{print $1}'`
   if test "${msg}" = "PASSED"; then
      # Compare succeed
      echo -e "Compare  : \033[32mPassed\033[0m ~> ${crun}"
   else
      # Compare failed
      echo -e "Compare  : \033[31mFailed\033[0m ~> ${crun}"
   fi
done

# Final time
final_time=`date --rfc-3339=seconds`
final_time_sec=`date +%s`
echo "Tests finished at ${final_time}" >> bump_ctest_log/execution.log

# Elapsed time
elapsed_time=$((final_time_sec-initial_time_sec))
echo "Elapsed time: ${elapsed_time} sec." >> bump_ctest_log/summary.log
echo "${stest} run tests succeed and ${ftest} failed" >> bump_ctest_log/summary.log
tail -n 2 bump_ctest_log/summary.log

exit 0
