#!/usr/bin/env bash
#----------------------------------------------------------------------
# Bash script: saber_doxygenation
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Directories
rootdir=$1/..
doc=${rootdir}/Documents
autodoc=${doc}/autodoc

# Languages
languages="*.cc *.f *.F90 *.h *.inc"

# Directory structure
dir=()
cd ${rootdir}
lev1=`ls -d */ 2> /dev/null`
for dir1 in ${lev1}; do
   cd ${dir1}
   lev2=`ls -d */ 2> /dev/null`
   for dir2 in ${lev2}; do
      cd ${dir2}
      list=`ls ${languages} 2> /dev/null`
      if test -n "${list}"; then
         dir+=("${dir1%?}/${dir2%?}")
      fi
      lev3=`ls -d */ 2> /dev/null`
      for dir3 in ${lev3}; do
         cd ${dir3}
         list=`ls ${languages} 2> /dev/null`
         if test -n "${list}"; then
            dir+=("${dir1%?}/${dir2%?}/${dir3%?}")
         fi
         lev4=`ls -d */ 2> /dev/null`
         for dir4 in ${lev4}; do
            cd ${dir4}
            list=`ls ${languages} 2> /dev/null`
            if test -n "${list}"; then
               dir+=("${dir1%?}/${dir2%?}/${dir3%?}/${dir4%?}")
            fi
            lev5=`ls -d */ 2> /dev/null`
            for dir5 in ${lev5}; do
               cd ${dir4}
               list=`ls ${languages} 2> /dev/null`
               if test -n "${list}"; then
                  dir+=("${dir1%?}/${dir2%?}/${dir3%?}/${dir4%?}/${dir5%?}")
               fi
            done
            cd ..
         done
         cd ..  
      done
      cd ..
   done
   cd ..
done

for index in ${!dir[*]}; do
   list=`ls ${rootdir}/${dir[$index]}/*.F90 ${rootdir}/${dir[$index]}/*.inc 2> /dev/null`

   for filename in ${list} ; do
      # Initialization
      new_type=false

      # While loop over lines
      echo ${filename}
      mv ${filename} ${filename}.old
      rm -f ${filename}
      OLDIFS=$IFS
      IFS=''
      cat ${filename}.old | while read -r line ; do
         done=false

         if [[ ${line} == "contains" ]] ; then
            new_type=false
         fi

         # Purpose
         if [[ ${line} == "! Purpose"* ]] ; then
            purpose=`echo ${line} | cut -c 12-`
            echo "!> "${purpose^} >> ${filename}
            done=true
         else
            # Intent
            if [[ ${line} == *"intent("* ]] ; then
               echo ${line} | sed -e s/"! "/"!< "/g  >> ${filename}
               done=true
            else
               # Type members
               if test "${new_type}" == "true" ; then
                  if [[ ${line} == *" :: "* ]] ; then
                     echo ${line} | sed -e s/"! "/"!< "/g  >> ${filename}
                     done=true
                  fi
               fi
            fi
         fi

         if [[ ${line} == "type "* ]] ; then
            new_type=true
         fi

         # Default
         if test "${done}" == "false" ; then
            echo ${line} >> ${filename}
         fi
      done
      IFS=$OLDIFS
      rm -f ${filename}.old
   done
done
