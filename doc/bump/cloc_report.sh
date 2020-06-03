#!/bin/bash
#----------------------------------------------------------------------
# Bash script: cloc_report
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Declare arrays
dir=(mains/bump mains/oops saber/bump saber/external saber/gaugrid saber/oops saber/util)
name=(mains_bump mains_oops saber_bump saber_external saber_gaugrid saber_oops saber_util)

if type "cloc" > /dev/null ; then
   # Directories
   src=$1/../../src
   doc=$1
   mkdir -p ${doc}

   # Get cloc report
   echo "Get cloc report"
   cd ..
   for index in ${!dir[*]}; do
      cloc --quiet --csv --exclude-lang=CMake --out=cloc_${name[$index]}.csv ${src}/${dir[$index]}
   done

   # Write reports
   echo "Refactor cloc report"
   echo -e "# cloc report\n" > ${doc}/cloc_report.md
   echo -e "Code report obtained with [CLOC](https://github.com/AlDanial/cloc).\n" >> ${doc}/cloc_report.md
   OLDIFS=$IFS
   IFS=,
   for index in ${!dir[*]}; do
      echo -e "## ${dir[$index]}\n" >> ${doc}/cloc_report.md
      i=0
      while read files language blank comment code dum ; do
         if test $i == 0 ; then
            ratio="$comment/$code ratio"
         else
            let ratio=100*comment/code
         fi
         echo -e "| $language | $files | $blank | $comment | $code | $ratio |" >> ${doc}/cloc_report.md
         if test $i == 0 ; then
            echo -e "|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|" >> ${doc}/cloc_report.md
         fi
         let i=i+1
      done < cloc_${name[$index]}.csv  
      echo -e "\n" >> ${doc}/cloc_report.md
   done
   IFS=$OLDIFS
   for index in ${!dir[*]}; do
      rm -f cloc_${name[$index]}.csv
   done
else
   echo "cloc not found: no cloc report"
fi
