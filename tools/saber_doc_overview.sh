#!/bin/bash
#----------------------------------------------------------------------
# Bash shell script: saber_doc_overview
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Directories
rootdir=$1/..
doc=${rootdir}/doc
autodoc=${doc}/autodoc

# Languages
languages="*.cc *.f *.F90 *.h"

# Title
echo "# SABER overview"
echo -e "# SABER overview" > ${doc}/overview.md

# Directory structure
echo "## Directory structure"
echo -e "## Directory structure" >> ${doc}/overview.md
echo -e "The SABER repository is organized as follows:" >> ${doc}/overview.md
dir=()
name=()
cd ${rootdir}
lev1=`ls -d */ 2> /dev/null`
for dir1 in ${lev1}; do
   cd ${dir1}
   desc=`cat .description 2> /dev/null`
   echo -e "* **"${dir1%?}"**: "${desc} >> ${doc}/overview.md
   list=`ls ${languages} 2> /dev/null`
   if test -n "${list}"; then
      dir+=("${dir1%?}")
      name+=("${dir1%?}")
   fi
   lev2=`ls -d */ 2> /dev/null`
   for dir2 in ${lev2}; do
      cd ${dir2}
      desc=`cat .description 2> /dev/null`
      echo -e "  * **"${dir2%?}"**: "${desc} >> ${doc}/overview.md
      list=`ls ${languages} 2> /dev/null`
      if test -n "${list}"; then
         dir+=("${dir1%?}/${dir2%?}")
         name+=("${dir1%?}_${dir2%?}")
      fi
      lev3=`ls -d */ 2> /dev/null`
      for dir3 in ${lev3}; do
         cd ${dir3}
         desc=`cat .description 2> /dev/null`
         echo -e "    * **"${dir3%?}"**: "${desc} >> ${doc}/overview.md
         list=`ls ${languages} 2> /dev/null`
         if test -n "${list}"; then
            dir+=("${dir1%?}/${dir2%?}/${dir3%?}")
            name+=("${dir1%?}_${dir2%?}_${dir3%?}")
         fi
         lev4=`ls -d */ 2> /dev/null`
         for dir4 in ${lev4}; do
            cd ${dir4}
            desc=`cat .description 2> /dev/null`
            echo -e "      * **"${dir4%?}"**: "${desc} >> ${doc}/overview.md
            list=`ls ${languages} 2> /dev/null`
            if test -n "${list}"; then
               dir+=("${dir1%?}/${dir2%?}/${dir3%?}/${dir4%?}")
               name+=("${dir1%?}_${dir2%?}_${dir3%?}_${dir4%?}")
            fi
            lev5=`ls -d */ 2> /dev/null`
            for dir5 in ${lev5}; do
               cd ${dir4}
               desc=`cat .description 2> /dev/null`
               echo -e "        * **"${dir5%?}"**: "${desc} >> ${doc}/overview.md
               list=`ls ${languages} 2> /dev/null`
               if test -n "${list}"; then
                  dir+=("${dir1%?}/${dir2%?}/${dir3%?}/${dir4%?}/${dir5%?}")
                  name+=("${dir1%?}_${dir2%?}_${dir3%?}_${dir4%?}_${dir5%?}")
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
echo -e "\n" >> ${doc}/overview.md

if type "cloc" > /dev/null ; then
   # Cloc report
   for index in ${!dir[*]}; do
      cloc --quiet --csv --exclude-lang=CMake --out=cloc_${name[$index]}.csv ${rootdir}/${dir[$index]}
   done

   # Code size and characteristics
   echo "## Code size and characteristics"
   echo -e "## Code size and characteristics" >> ${doc}/overview.md
   echo -e "Code report obtained with [CLOC](https://github.com/AlDanial/cloc).\n" >> ${doc}/overview.md
   OLDIFS=$IFS
   IFS=,
   for index in ${!dir[*]}; do
      echo -e "### ${dir[$index]}\n" >> ${doc}/overview.md
      i=0
      while read files language blank comment code dum ; do
         if test $i == 0 ; then
            ratio="${comment}/${code} ratio"
         else
            let ratio=100*comment/code
            ratio="${ratio} %"
         fi
         echo -e "| ${language} | ${files} | ${blank} | ${comment} | ${code} | ${ratio} |" >> ${doc}/overview.md
         if test $i == 0 ; then
            echo -e "|:--------:|:--------:|:--------:|:--------:|:--------:|:--------:|" >> ${doc}/overview.md
         fi
         let i=i+1
      done < cloc_${name[$index]}.csv  
      echo -e "\n" >> ${doc}/overview.md
   done
   IFS=$OLDIFS
   for index in ${!dir[*]}; do
      rm -f cloc_${name[$index]}.csv
   done
else
   echo "cloc not found: no cloc report"
fi

# Code auto-documentation
echo "## Code auto-documentation"
echo -e "## Code auto-documentation" >> ${doc}/overview.md
for index in ${!dir[*]}; do
   list=`ls ${rootdir}/${dir[$index]}/*.F90 2> /dev/null`
   desc=`cat ${rootdir}/${dir[$index]}/.description 2> /dev/null`

   echo -e "### ${dir[$index]}: ${desc}\n" >> ${doc}/overview.md
   echo -e "| Name | Purpose |" >> ${doc}/overview.md
   echo -e "| :--: | :---------- |" >> ${doc}/overview.md

   for filename in ${list} ; do
      # Initialization
      new_module=false
      new_purpose=false
      new_subfunc=false
      type_bound=false
      i=-2

      # While loop over lines
      while IFS= read -r line ; do
         # Get keywords
         word=`echo ${line} | cut -c -10`
         if test "${word}" = "! Module: " ; then
            module=`echo ${line} | cut -c 11-`
            if test "${category}" = "derived_type" ; then
               class=${module#*type_}
            fi
            new_module=true
         fi
         if test "${word}" = "! Purpose:" ; then
            purpose=`echo ${line} | cut -c 12-`
            new_purpose=true
         fi
         if test "${word}" = "! Subrouti" ; then
            subfunc=`echo ${line} | cut -c 15-`
            if test "${category}" = "derived_type" ; then
               if test "${subfunc#$class*}" != "${subfunc}" ; then
                  type_bound=true
               fi
            fi
            subfunc_type=subroutine
            new_subfunc=true
         fi
         if test "${word}" = "! Function" ; then
            subfunc=`echo ${line} | cut -c 13-`
            if test "${category}" = "derived_type" ; then
               if test "${subfunc#${class}*}" != "${subfunc}" ; then
                  type_bound=true
               fi
            fi
            subfunc_type=function
            new_subfunc=true
         fi

         # Increment line index
         let i=i+1

         if test "${new_purpose}" = "true" ; then
            # New module
            if test "${new_module}" = "true" ; then
               echo -e "| [${module}](autodoc/${module}.md) | ${purpose} |" >> ${doc}/overview.md

               type_bound=false
               cat<<EOFMOD > ${autodoc}/${module}.md
# Module ${module}

| Type | Name | Purpose |
| :--: | :--: | :---------- |
EOFMOD
               new_module=false
            fi

            # New subroutine/function
            if test "${new_subfunc}" = "true" ; then
               if test "${type_bound}" = "true" ; then
                  echo -e "| ${subfunc_type} | [${class}\%] [${subfunc#${class}_}](https://github.com/JCSDA/saber/tree/develop/${dir[$index]}/${module}.F90#L${i}) | ${purpose} |" >> ${autodoc}/${module}.md
               else
                  echo -e "| ${subfunc_type} | [${subfunc}](https://github.com/JCSDA/saber/tree/develop/${dir[$index]}/${module}.F90#L${i}) | ${purpose} |" >> ${autodoc}/${module}.md
               fi
               new_subfunc=false
               type_bound=false
            fi
         fi

         # Reset
         new_purpose=false
      done < ${filename}
   done
   echo -e "\n" >> ${doc}/overview.md
done
