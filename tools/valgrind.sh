#!/bin/sh
#----------------------------------------------------------------------
# Shell script: valgrind
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Parameters
exe=$1
input=$2
output=$3

# Initialize exit status
status=0

if type "valgrind --help" > /dev/null ; then
   # Run valgrind
   valgrind --log-file="${output}/valgrind.out" -q ${exe} ${input} ${output}

   # Loop over file prefixes
   for prefix in "bump_" "type_" "tools_" ; do
      # Get error output, find files starting with a given prefix, select column with file and line, remove duplicates
      errorlist=`grep "  at" ${output}/valgrind.out | grep __$prefix | awk '{print $5}' | awk '{for (i=1;i<=NF;i++) if (!a[$i]++) printf("%s%s",$i,FS)}{printf("\n")}'`
      for error in $errorlist ; do
         # Remove parentheses
         error=${error#(*}
         error=${error%)}

         # Get file and line
         file=${error%:*}
         line=${error#*:}

         # Print results
         echo "\e[31mError in file ${file} at line ${line}\e[0m"
         status=2
      done
   done
else
   echo "\e[31mCannot find command: valgrind\e[0m"
   status=2
fi

# Exit
exit ${status}
