#!/usr/bin/env bash
#----------------------------------------------------------------------
# Bash script: saber_parallel
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

function ProgressBar {
   let _progress=(${1}*100/${2}*100)/100
   let _done=(${_progress}*4)/10
   let _left=40-$_done
   _fill=$(printf "%${_done}s")
   _empty=$(printf "%${_left}s")
   string=`printf "Progress : [${_fill// /#}${_empty// /-}] ${_progress}%% (${3})"`
   echo -e "\e[1A\e[K${string}"
}

function PrintFailed {
   echo -e "\e[1A\e[K\033[31m${1} failed\033[0m"
   echo -e ""
}

# Check bash version
if test ${BASH_VERSINFO:-0} -lt 4; then
   echo "Bash version lower than 4 does not support associative arrays"
   exit 1
fi

# Check location
PWD=`pwd`
if test "${PWD##*/}" != "test"; then
   echo "This script should be run from \${build_directory}/saber/test"
   exit 1
fi
tmp=${PWD%/test}
if test "${tmp##*/}" != "saber"; then
   echo "This script should be run from \${build_directory}/saber/test"
   exit 1
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
   exit 1
fi
if test ${nproc} -gt ${nprocmax}; then
   echo "Only ${nprocmax} processor(s) are available"
   exit 1
fi
echo "Tests run on ${nproc} logical core(s) (${nprocmax} available)"

# Compiler
mpirun_path=`which mpirun`
if [[ ${mpirun_path} == *"openmpi"* ]]; then
   list_command="-cpu-list "
fi
if [[ ${mpirun_path} == *"intel"* ]]; then
   list_command="-genv I_MPI_PIN_PROCESSOR_LIST="   
fi
if test "${list_command}" = "" ; then
   echo "Cannot find what compiler is used for mpirun"
   exit 1
fi

# Make temporary directory
mkdir -p saber_ctest_log
rm -fr saber_ctest_log/*

# Initial time
initial_time=`date`
initial_time_sec=`date +%s`
echo "Tests start at ${initial_time}" > saber_ctest_log/execution.log

# List tests
list_get=`ctest -N | grep ": get_" | awk '{print $(NF)}'`
list_run=`ctest -N | grep ": test_bump" | grep _run | awk '{print $(NF)}'`
list_compare=`ctest -N | grep ": test_bump" | grep _compare | awk '{print $(NF)}'`
list_post=`ctest -N | grep ": test_bump" | grep _post | awk '{print $(NF)}'`
list_plot=`ctest -N | grep ": test_bump" | grep _plot | awk '{print $(NF)}'`
list_valgrind=`ctest -N | grep ": test_bump" | grep _valgrind | awk '{print $(NF)}'`
list_qg=`ctest -N | grep ": test_qg" | awk '{print $(NF)}'`
list_interpolation=`ctest -N | grep ": test_interpolation" | awk '{print $(NF)}'`

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
list_post_array=(${list_post})
ntest_post=${#list_post_array[@]}
stest_post=0
ftest_post=0
list_plot_array=(${list_plot})
ntest_plot=${#list_plot_array[@]}
stest_plot=0
ftest_plot=0
list_valgrind_array=(${list_valgrind})
ntest_valgrind=${#list_valgrind_array[@]}
stest_valgrind=0
ftest_valgrind=0
list_qg_array=(${list_qg})
ntest_qg=${#list_qg_array[@]}
stest_qg=0
ftest_qg=0
list_interpolation_array=(${list_interpolation})
ntest_interpolation=${#list_interpolation_array[@]}
stest_interpolation=0
ftest_interpolation=0
ntest=$((ntest_run+ntest_compare+ntest_post+ntest_plot+ntest_valgrind+ntest_qg+ntest_interpolation))
itest=0
if test ${ntest} = 0; then
   echo "No test detected, this script should be run from \${build_directory}/saber/test"
   exit 1
fi

# Declare pids array
declare -A pids=()

# Download tests
for get in ${list_get}; do
   echo "Handling process ${get}" >> saber_ctest_log/execution.log
   ctest -R ${get}\$ > saber_ctest_log/${get}.log 2> saber_ctest_log/${get}.err &
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
echo -e ""
ProgressBar 0 ${ntest} ""

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
   ctest --timeout 0.00001 -VV -R ${run}\$ > saber_ctest_log/${run}_fail.log 2>/dev/null
   tmp=`grep "Test command:" saber_ctest_log/${run}_fail.log`
   exe=`echo ${tmp} | awk '{print $7,$8,$9}'`
   tmp1=${run%%_run}
   omp=${tmp1##*-}
   tmp2=${tmp1%-*}
   mpi=${tmp2##*_}
   mpixomp=$((mpi*omp))
   if [ ${mpixomp} -gt ${nproc} ]; then
      echo "Too many CPUs are required for process ${run} (${mpixomp} > ${nproc})"
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
            if test "${err}" = "0"; then
               # Run passed
               echo "${crun} passed" >> saber_ctest_log/execution.log
               stest_run=$((stest_run+1))
            else
               # Run failed
               echo "${crun} failed" >> saber_ctest_log/execution.log
               ftest_run=$((ftest_run+1))
               PrintFailed ${crun}
            fi
            ProgressBar ${itest} ${ntest} ${crun}
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
   full_cmd="mpirun ${list_command}${cpu_list} -n ${mpi} ${exe} > saber_ctest_log/${run}.log 2> saber_ctest_log/${run}.err &"
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
         percent=$((100*${itest}/${ntest}))

         if test "${err}" = "0"; then
            # Run passed
            echo "${crun} passed" >> saber_ctest_log/execution.log
            stest_run=$((stest_run+1))
         else
            # Run failed
            echo "${crun} failed" >> saber_ctest_log/execution.log
            ftest_run=$((ftest_run+1))
            PrintFailed ${crun}
         fi
         ProgressBar ${itest} ${ntest} ${crun}
      fi
   done
done

# Re-initialize cpus array
for i in $(seq 1 ${nproc}); do
   cpus[$i]=""
done
navail=${nproc}

# Compare tests
for compare in ${list_compare}; do
   echo "Handling process ${compare}" >> saber_ctest_log/execution.log

   # Get command and arguments
   ctest --timeout 0.00001 -VV -R ${compare}\$ > saber_ctest_log/${compare}_fail.log 2>/dev/null
   tmp=`grep "Test command:" saber_ctest_log/${compare}_fail.log`
   exe=`echo ${tmp} | awk '{print $4,$5,$6,$7,$8}'`
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
            if test "${err}" = "0"; then
               # Compare passed
               echo "${ccompare} passed" >> saber_ctest_log/execution.log
               stest_compare=$((stest_compare+1))
            else
               # Compare failed
               echo "${ccompare} failed" >> saber_ctest_log/execution.log
               ftest_compare=$((ftest_compare+1))
               PrintFailed ${ccompare}
            fi
            ProgressBar ${itest} ${ntest} ${ccompare}
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
   full_cmd="mpirun ${list_command}\"${cpu_list}\" \"-n\" \"${mpi}\" ${exe} > saber_ctest_log/${compare}.log 2> saber_ctest_log/${compare}.err &"
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
         if test "${err}" = "0"; then
            # Compare passed
            echo "${ccompare} passed" >> saber_ctest_log/execution.log
            stest_compare=$((stest_compare+1))
         else
            # Compare failed
            echo "${ccompare} failed" >> saber_ctest_log/execution.log
            ftest_compare=$((ftest_compare+1))
            PrintFailed ${ccompare}
         fi
         ProgressBar ${itest} ${ntest} ${ccompare}
      fi
   done
done

# Re-initialize cpus array
for i in $(seq 1 ${nproc}); do
   cpus[$i]=""
done
navail=${nproc}

# Post tests
for post in ${list_post}; do
   echo "Handling process ${post}" >> saber_ctest_log/execution.log

   # Get command and arguments
   ctest --timeout 0.00001 -VV -R ${post}\$ > saber_ctest_log/${post}_fail.log 2>/dev/null
   tmp=`grep "Test command:" saber_ctest_log/${post}_fail.log`
   exe=`echo ${tmp} | awk '{print $4,$5,$6,$7,$8}'`
   mpi=1
   omp=1
   mpixomp=$((mpi*omp))
   rm -f saber_ctest_log/${post}_fail.log

   # Print info
   echo "   At this time, ${navail} CPUs are available and ${mpixomp} are required" >> saber_ctest_log/execution.log

   # Wait for available CPUs
   while [ ${navail} -lt ${mpixomp} ]; do
      # Sleep
      sleep 0.01

      # Loop over running processes
      for cpost in "${!pids[@]}"; do
         # Check if the process is not running
         if ! kill -0 ${pids[${cpost}]} 2>/dev/null; then
            # Unset process from the list
            unset pids[${cpost}]

            # CPUs for this process are now available
            navail=0
            for i in "${!cpus[@]}"; do
               if test "${cpus[${i}]}" = ""; then
                  navail=$((navail+1))
               fi
               if test "${cpus[${i}]}" = "${cpost}"; then
                  # This CPU is now available
                  cpus[${i}]=""
                  navail=$((navail+1))
               fi
            done
            echo "   Process ${cpost} is done, ${navail} are now are available" >> saber_ctest_log/execution.log

            # Check if this process passed
            err=`wc -l saber_ctest_log/${cpost}.err | awk '{print $1}'`
            itest=$((itest+1))
            if test "${err}" = "0"; then
               # Post passed
               echo "${cpost} passed" >> saber_ctest_log/execution.log
               stest_post=$((stest_post+1))
            else
               # Post failed
               echo "${cpost} failed" >> saber_ctest_log/execution.log
               ftest_post=$((ftest_post+1))
               PrintFailed ${cpost}
            fi
            ProgressBar ${itest} ${ntest} ${cpost}
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
            cpus[${i}]="${post}"
         fi
      fi
   done

   # Fire!
   echo "   Process ${post} is launched on CPUs ${cpu_list}" >> saber_ctest_log/execution.log
   export OMP_NUM_THREADS=${omp}
   full_cmd="mpirun ${list_command}\"${cpu_list}\" \"-n\" \"${mpi}\" ${exe} > saber_ctest_log/${post}.log 2> saber_ctest_log/${post}.err &"
   eval ${full_cmd}
   pids[${post}]=$!
done

# Wait for all processes to finish
while [ ${#pids[@]} -gt 0 ]; do
   # Sleep
   sleep 0.01

   # Loop over running processes
   for cpost in "${!pids[@]}"; do
      # Check if the process is not running
      if ! kill -0 ${pids[${cpost}]} 2>/dev/null; then
         # Unset process from the list
         unset pids[${cpost}]

         # Check if this process passed
         err=`wc -l saber_ctest_log/${cpost}.err | awk '{print $1}'`
         itest=$((itest+1))
         if test "${err}" = "0"; then
            # Post passed
            echo "${cpost} passed" >> saber_ctest_log/execution.log
            stest_post=$((stest_post+1))
         else
            # Post failed
            echo "${cpost} failed" >> saber_ctest_log/execution.log
            ftest_post=$((ftest_post+1))
            PrintFailed ${cpost}
         fi
         ProgressBar ${itest} ${ntest} ${cpost}
      fi
   done
done

# Re-initialize cpus array
for i in $(seq 1 ${nproc}); do
   cpus[$i]=""
done
navail=${nproc}

# Plot tests
for plot in ${list_plot}; do
   echo "Handling process ${plot}" >> saber_ctest_log/execution.log

   # Get command and arguments
   ctest --timeout 0.00001 -VV -R ${plot}\$ > saber_ctest_log/${plot}_fail.log 2>/dev/null
   tmp=`grep "Test command:" saber_ctest_log/${plot}_fail.log`
   exe=`echo ${tmp} | awk '{print $4,$5,$6,$7,$8,$9}'`
   mpi=1
   omp=1
   mpixomp=$((mpi*omp))
   rm -f saber_ctest_log/${plot}_fail.log

   # Print info
   echo "   At this time, ${navail} CPUs are available and ${mpixomp} are required" >> saber_ctest_log/execution.log

   # Wait for available CPUs
   while [ ${navail} -lt ${mpixomp} ]; do
      # Sleep
      sleep 0.01

      # Loop over running processes
      for cplot in "${!pids[@]}"; do
         # Check if the process is not running
         if ! kill -0 ${pids[${cplot}]} 2>/dev/null; then
            # Unset process from the list
            unset pids[${cplot}]

            # CPUs for this process are now available
            navail=0
            for i in "${!cpus[@]}"; do
               if test "${cpus[${i}]}" = ""; then
                  navail=$((navail+1))
               fi
               if test "${cpus[${i}]}" = "${cplot}"; then
                  # This CPU is now available
                  cpus[${i}]=""
                  navail=$((navail+1))
               fi
            done
            echo "   Process ${cplot} is done, ${navail} are now are available" >> saber_ctest_log/execution.log

            # Check if this process passed
            err=`wc -l saber_ctest_log/${cplot}.err | awk '{print $1}'`
            itest=$((itest+1))
            if test "${err}" = "0"; then
               # Plot passed
               echo "${cplot} passed" >> saber_ctest_log/execution.log
               stest_plot=$((stest_plot+1))
            else
               # Plot failed
               echo "${cplot} failed" >> saber_ctest_log/execution.log
               ftest_plot=$((ftest_plot+1))
               PrintFailed ${cplot}
            fi
            ProgressBar ${itest} ${ntest} ${cplot}
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
            cpus[${i}]="${plot}"
         fi
      fi
   done

   # Fire!
   echo "   Process ${plot} is launched on CPUs ${cpu_list}" >> saber_ctest_log/execution.log
   export OMP_NUM_THREADS=${omp}
   full_cmd="mpirun \"-n\" \"${mpi}\" \"-cpu-list\" \"${cpu_list}\" ${exe} > saber_ctest_log/${plot}.log 2> saber_ctest_log/${plot}.err &"
   eval ${full_cmd}
   pids[${plot}]=$!
done

# Wait for all processes to finish
while [ ${#pids[@]} -gt 0 ]; do
   # Sleep
   sleep 0.01

   # Loop over running processes
   for cplot in "${!pids[@]}"; do
      # Check if the process is not running
      if ! kill -0 ${pids[${cplot}]} 2>/dev/null; then
         # Unset process from the list
         unset pids[${cplot}]

         # Check if this process passed
         err=`wc -l saber_ctest_log/${cplot}.err | awk '{print $1}'`
         itest=$((itest+1))
         if test "${err}" = "0"; then
            # Plot passed
            echo "${cplot} passed" >> saber_ctest_log/execution.log
            stest_plot=$((stest_plot+1))
         else
            # Plot failed
            echo "${cplot} failed" >> saber_ctest_log/execution.log
            ftest_plot=$((ftest_plot+1))
            PrintFailed ${cplot}
         fi
         ProgressBar ${itest} ${ntest} ${cplot}
      fi
   done
done

# Re-initialize cpus array
for i in $(seq 1 ${nproc}); do
   cpus[$i]=""
done
navail=${nproc}

# Valgrind tests
for valgrind in ${list_valgrind}; do
   echo "Handling process ${valgrind}" >> saber_ctest_log/execution.log

   # Get command and arguments
   ctest --timeout 0.00001 -VV -R ${valgrind}\$ > saber_ctest_log/${valgrind}_fail.log 2>/dev/null
   tmp=`grep "Test command:" saber_ctest_log/${valgrind}_fail.log`
   cmd=`echo ${tmp} | awk '{print $4}'`
   exe=`echo ${tmp} | awk '{print $7,$8,$9,$10}'`
   omp=1
   mpi=1
   mpixomp=$((mpi*omp))
   rm -f saber_ctest_log/${valgrind}_fail.log

   # Print info
   echo "   At this time, ${navail} CPUs are available and ${mpixomp} are required" >> saber_ctest_log/execution.log

   # Wait for available CPUs
   while [ ${navail} -lt ${mpixomp} ]; do
      # Sleep
      sleep 0.01

      # Loop over running processes
      for cvalgrind in "${!pids[@]}"; do
         # Check if the process is not running
         if ! kill -0 ${pids[${cvalgrind}]} 2>/dev/null; then
            # Unset process from the list
            unset pids[${cvalgrind}]

            # CPUs for this process are now available
            navail=0
            for i in "${!cpus[@]}"; do
               if test "${cpus[${i}]}" = ""; then
                  navail=$((navail+1))
               fi
               if test "${cpus[${i}]}" = "${cvalgrind}"; then
                  # This CPU is now available
                  cpus[${i}]=""
                  navail=$((navail+1))
               fi
            done
            echo "   Process ${cvalgrind} is done, ${navail} are now are available" >> saber_ctest_log/execution.log

            # Check if this process passed
            err=`wc -l saber_ctest_log/${cvalgrind}.err | awk '{print $1}'`
            itest=$((itest+1))
            if test "${err}" = "0"; then
               # Run passed
               echo "${cvalgrind} passed" >> saber_ctest_log/execution.log
               stest_valgrind=$((stest_valgrind+1))
            else
               # Run failed
               echo "${cvalgrind} failed" >> saber_ctest_log/execution.log
               ftest_valgrind=$((ftest_valgrind+1))
               PrintFailed ${cvalgrind}
            fi
            ProgressBar ${itest} ${ntest} ${cvalgrind}
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
            cpus[${i}]="${valgrind}"
         fi
      fi
   done

   # Fire!
   echo "   Process ${valgrind} is launched on CPUs ${cpu_list}" >> saber_ctest_log/execution.log
   export OMP_NUM_THREADS=${omp}
   full_cmd="${cmd} \"-n\" \"${mpi}\" \"-cpu-list\" \"${cpu_list}\" ${exe} > saber_ctest_log/${valgrind}.log 2> saber_ctest_log/${valgrind}.err &"
   eval ${full_cmd}
   pids[${valgrind}]=$!
done

# Wait for all processes to finish
while [ ${#pids[@]} -gt 0 ]; do
   # Sleep
   sleep 0.01

   # Loop over running processes
   for cvalgrind in "${!pids[@]}"; do
      # Check if the process is not running
      if ! kill -0 ${pids[${cvalgrind}]} 2>/dev/null; then
         # Unset process from the list
         unset pids[${cvalgrind}]

         # Check if this process passed
         err=`wc -l saber_ctest_log/${cvalgrind}.err | awk '{print $1}'`
         itest=$((itest+1))
         if test "${err}" = "0"; then
            # Run passed
            echo "${cvalgrind} passed" >> saber_ctest_log/execution.log
            stest_valgrind=$((stest_valgrind+1))
         else
            # Run failed
            echo "${cvalgrind} failed" >> saber_ctest_log/execution.log
            ftest_valgrind=$((ftest_valgrind+1))
            PrintFailed ${cvalgrind}
         fi
         ProgressBar ${itest} ${ntest} ${cvalgrind}
      fi
   done
done

# Re-initialize cpus array
for i in $(seq 1 ${nproc}); do
   cpus[$i]=""
done
navail=${nproc}

# QG tests
for qg in ${list_qg}; do
   echo "Handling process ${qg}" >> saber_ctest_log/execution.log

   # Get command and arguments
   ctest -VV -R ${qg}\$ > saber_ctest_log/${qg}.log 2> saber_ctest_log/${qg}.err

   # Check if this process passed
   err=`wc -l saber_ctest_log/${qg}.err | awk '{print $1}'`
   itest=$((itest+1))
   if test "${err}" = "0"; then
      # QG passed
      echo "${qg} passed" >> saber_ctest_log/execution.log
      stest_qg=$((stest_qg+1))
   else
      # QG failed
      echo "${qg} failed" >> saber_ctest_log/execution.log
      ftest_qg=$((ftest_qg+1))
      PrintFailed ${qg}
   fi
   ProgressBar ${itest} ${ntest} ${qg}
done

# Interpolation tests
for interpolation in ${list_interpolation}; do
   echo "Handling process ${interpolation}" >> saber_ctest_log/execution.log

   # Get command and arguments
   ctest -VV -R ${interpolation}\$ > saber_ctest_log/${interpolation}.log 2> saber_ctest_log/${interpolation}.err

   # Check if this process passed
   err=`wc -l saber_ctest_log/${interpolation}.err | awk '{print $1}'`
   itest=$((itest+1))
   if test "${err}" = "0"; then
      # Interpolation passed
      echo "${interpolation} passed" >> saber_ctest_log/execution.log
      stest_interpolation=$((stest_interpolation+1))
   else
      # Interpolation failed
      echo "${interpolation} failed" >> saber_ctest_log/execution.log
      ftest_interpolation=$((ftest_interpolation+1))
      PrintFailed ${interpolation}
   fi
   ProgressBar ${itest} ${ntest} ${interpolation}
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

# Print summary
echo -e ""
echo -e "Summary:"

if test "${ntest_get}" -gt "0"; then
   # Download tests
   stest=`printf "%03d" $((stest_get))`
   ftest=`printf "%03d" $((ftest_get))`
   echo -e "  Download:      \033[32m${stest}\033[0m tests passed and \033[31m${ftest}\033[0m failed"
fi

if test "${ntest_run}" -gt "0"; then
   # Run tests
   stest=`printf "%03d" $((stest_run))`
   ftest=`printf "%03d" $((ftest_run))`
   echo -e "  Run:           \033[32m${stest}\033[0m tests passed and \033[31m${ftest}\033[0m failed"
fi

if test "${ntest_compare}" -gt "0"; then
   # Compare tests
   stest=`printf "%03d" $((stest_compare))`
   ftest=`printf "%03d" $((ftest_compare))`
   echo -e "  Compare:       \033[32m${stest}\033[0m tests passed and \033[31m${ftest}\033[0m failed"
fi

if test "${ntest_post}" -gt "0"; then
   # post tests
   stest=`printf "%03d" $((stest_post))`
   ftest=`printf "%03d" $((ftest_post))`
   echo -e "  Post:          \033[32m${stest}\033[0m tests passed and \033[31m${ftest}\033[0m failed"
fi

if test "${ntest_plot}" -gt "0"; then
   # Plot tests
   stest=`printf "%03d" $((stest_plot))`
   ftest=`printf "%03d" $((ftest_plot))`
   echo -e "  Plot:          \033[32m${stest}\033[0m tests passed and \033[31m${ftest}\033[0m failed"
fi

if test "${ntest_valgrind}" -gt "0"; then
   # Valgrind tests
   stest=`printf "%03d" $((stest_valgrind))`
   ftest=`printf "%03d" $((ftest_valgrind))`
   echo -e "  Valgrind:      \033[32m${stest}\033[0m tests passed and \033[31m${ftest}\033[0m failed"
fi

if test "${ntest_qg}" -gt "0"; then
   # QG tests
   stest=`printf "%03d" $((stest_qg))`
   ftest=`printf "%03d" $((ftest_qg))`
   echo -e "  QG:            \033[32m${stest}\033[0m tests passed and \033[31m${ftest}\033[0m failed"
fi

if test "${ntest_interpolation}" -gt "0"; then
   # Interpolation tests
   stest=`printf "%03d" $((stest_interpolation))`
   ftest=`printf "%03d" $((ftest_interpolation))`
   echo -e "  Interpolation: \033[32m${stest}\033[0m tests passed and \033[31m${ftest}\033[0m failed"
fi

exit 0
