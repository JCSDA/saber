#!/bin/bash
# This bash script can be used to evaluate the code coverage offline (not within CI testing).
# First, you need to set the environment variable ENABLE_OFFLINE_CODECOV at ON at build time,
# in order to activate the use of gcov and compute the code coverage files.
# Second, compile the code (in debug mode).
# Third, run the tests using ctest.
# Finally, update this script with your own build directory ('build_dir' variable) and run it to
# generate an HTML summary of the results (requires lcov: https://github.com/linux-test-project/lcov).
# The optional argument 'final_tar' can be used to produce a tar file with the final results instead of
# viewing them with firefox.
#
# Usage: saber_gcov.sh [final_tar]

if [[ $# -gt 1 ]]; then
   echo "Error: wrong number of arguments"
   echo "Usage: saber_gcov.sh [final_tar]"
   echo "- final_tar is optional and can be true or false (default is false), if set to true, the final results are gathered into a tar file"
   exit 1
fi

# Input paramters
if [[ $# -eq 1 ]]; then
   final_tar=$1
else
   final_tar=false
fi

# Hard-coded parameters
build_dir=${HOME}/build/gnu_10.3.0/bundle_debug
saber_dirs="
bump
external
gsi
interpolation
oops
spectralb
util
vader"

# Process output with gcov
for dir in ${saber_dirs}; do   
   cd ${build_dir}/saber/src/saber/CMakeFiles/saber.dir/${dir}
   for file in $(ls *.F90.gcda *.F90.gcno); do
      newfile=$(echo ${file} | sed -e 's/\.F90\././1')
      if [ ! -f ${newfile} ]; then
        ln -sf ${file} ${newfile}
      fi
   done
   ln -sf ${build_dir}/saber/src/saber/${dir}/*.F90 .
   gcov *.F90
done

# Process output with lcov
cd ${build_dir}/saber/src/saber/CMakeFiles/saber.dir
rm -fr coverage.info html
lcov --gcov-tool gcov --capture --directory . --output-file coverage.info

# Generate html
genhtml --output-directory html coverage.info

if [[ ${final_tar} == "true" ]]; then
   # Tar results
   tar -cf coverage.tar html
elif [[ ${final_tar} == "false" ]]; then
   # Open html with firefox
   firefox html/saber/src/saber/bump/index.html
else
   # Wrong final_tar argument
   echo "Error: wrong final_tar argument, should be true or false"
   exit 2
fi
