#!/bin/ksh
#----------------------------------------------------------------------
# Korn shell script: cloc_report
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

if type "cloc" > /dev/null ; then
   # Directories
   src=$1/../../src/bump
   doc=$1
   mkdir -p ${doc}

   # Get cloc report
   echo "       Get cloc report"
   cd ..
   cloc --exclude-dir=external --quiet --csv --out=cloc.csv ${src}
   cloc --quiet --csv --out=cloc_external.csv ${src}/external
   
   # Write report
   echo "       Refactor cloc report"
   printf "# CLOC_REPORT\n\n" > ${doc}/CLOC_REPORT.md
   printf "Code report obtained with [CLOC](https://github.com/AlDanial/cloc).\n\n" >> ${doc}/CLOC_REPORT.md
   OLDIFS=$IFS
   IFS=,
   printf "**BUMP code:** \n\n" >> ${doc}/CLOC_REPORT.md
   i=0
   while read files language blank comment code dum ; do
      if test $i == 0 ; then
         ratio="$comment/$code ratio"
      else
         let ratio=100*comment/code
      fi
      printf "| $language | $files | $blank | $comment | $code | $ratio |\n" >> ${doc}/CLOC_REPORT.md
      if test $i == 0 ; then
         printf "|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|\n" >> ${doc}/CLOC_REPORT.md
      fi
      let i=i+1
   done < cloc.csv
   printf "\n" >> ${doc}/CLOC_REPORT.md
   printf "**External code:** \n\n" >> ${doc}/CLOC_REPORT.md
   i=0
   while read files language blank comment code dum ; do
      if test $i == 0 ; then
         ratio="$comment/$code ratio"
      else
         let ratio=100*comment/code
      fi
      printf "| $language | $files | $blank | $comment | $code | $ratio |\n" >> ${doc}/CLOC_REPORT.md
      if test $i == 0 ; then
         printf "|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|\n" >> ${doc}/CLOC_REPORT.md
      fi
      let i=i+1
   done < cloc_external.csv
   IFS=$OLDIFS
   rm -f cloc.csv cloc_external.csv
else
   echo "cloc not found: no cloc report"
fi
