#!/bin/bash
# This script can be used to update reference files in SABER. Update this script with your own
# source and build directories ('src_dir' and 'build_dir' variables, respectively), and run it.
# The 'rerun' variable can be set to 1 to rerun tests after the reference update and make sure
# they are actually passing. Optional argument are the indices of the first and last tests to
# execute (test indices can be found with "ctest -N").
#
# Usage: saber_update.sh [first_test_index] [last_test_index]

# Hard-coded parameters
src_dir=${HOME}/code/bundle/saber
build_dir=${HOME}/build/gnu_10.3.0/bundle_RelWithDebInfo
rerun=1

# Go to SABER build directory
cd ${build_dir}/saber

indices=`ctest -N | grep saber_test_ | gawk '{print $2}'`
for index in ${indices}; do
   name=`ctest -N | grep "${index}" | gawk '{print $3}'`
   index=${index##\#}
   index=${index%%:*}

   if [[ $# -eq 0 ]]; then
      valid=true
   elif [[ $# -eq 1 ]]; then
      if [[ ${index} -ge $1 ]]; then
         valid=true
      else
         valid=false
      fi
   elif [[ $# -eq 2 ]]; then
      if [[ ${index} -ge $1 ]] && [[ ${index} -le $2 ]]; then
         valid=true
      else
         valid=false
      fi
   else
      echo "Error: wrong number of arguments"
      echo "Usage: saber_update.sh [first_test_index] [last_test_index]"
      exit 1
   fi

   if [[ ${valid} == "true" ]]; then
      echo "Test #"${index}": "${name}
      ctest -VV -I ${index},${index} > tmp.log.out 2>/dev/null
      passed=`grep -s "100% tests passed" tmp.log.out | wc -l`
      if [[ ${passed} -eq 0 ]]; then
         isref=`grep -s "Comparing to reference file:" tmp.log.out | wc -l`
         if [[ ${isref} -eq "1" ]]; then
            ref=`grep "Comparing to reference file:" tmp.log.out | gawk '{print $7}'`
            echo "  Reference file: "${src_dir}/test/${ref}
            length=`grep -si "Test     :" tmp.log.out | gawk '{print $1}' | uniq -c | gawk '{print $2}' | gawk '{print length}'`
            rem=$((13+length))
            grep -si "Test     :" tmp.log.out | cut -c ${rem}- > ${src_dir}/test/${ref}
         else
            echo "Cannot find reference file"
            cat tmp.log.out
            exit 1
         fi
         if [[ ${rerun} -eq 1 ]]; then
            ctest -I ${index},${index} > tmp.log.out 2>/dev/null
            success=`grep -s "Passed" tmp.log.out | wc -l`
            if [[ ${success} -eq "1" ]] ; then
               echo -e "  Rerun \e[32mpassed\e[0m"
            else
               echo -e "  Rerun \e[31mfailed\e[0m"
            fi
         fi
      else
         echo -e "  Run \e[32mpassed\e[0m"
      fi
      rm -f tmp.log.out
   fi
done

exit 0
