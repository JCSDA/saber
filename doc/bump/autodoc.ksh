#!/bin/ksh
#----------------------------------------------------------------------
# Korn shell script: autodoc
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Directories
src=$1/../../src/saber/bump
doc=$1
mkdir -p ${doc}
autodoc=$1/autodoc
mkdir -p ${autodoc}

# Introduction
cat<<EOFINTRO > ${doc}/code_autodoc.md
# Code auto-documentation

The source code is organized in three categories:

1. Useful tools for the whole code: tools_[...].F90
2. Derived types and associated methods: type_[...].F90
3. External tools: external/[...].F90

| Category | Name | Purpose |
| :------: | :--: | :---------- |
EOFINTRO


for category in "tools" "derived_type" "external_tools" ; do
   echo "       Category: "${category}
   if test "${category}" = "main" ; then
      list=${mains}/bump.F90

   fi
   if test "${category}" = "tools" ; then
      list=`ls ${src}/tools_*.F90`
   fi
   if test "${category}" = "derived_type" ; then
      list=`ls ${src}/type_*.F90`
   fi
   if test "${category}" = "external_tools" ; then
      list=`ls ${src}/external/*.F90`
   fi

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
            echo "          Module: "${module}
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
               printf "| ${category} | [${module}](autodoc/${module}.md) | ${purpose} |\n" >> ${doc}/code_autodoc.md

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
                  printf "| ${subfunc_type} | [${class}\%] [${subfunc#${class}_}](https://github.com/JCSDA/saber/src/saber/bump/${module}.F90#L${i}) | ${purpose} |\n" >> ${autodoc}/${module}.md
               else
                  printf "| ${subfunc_type} | [${subfunc}](https://github.com/JCSDA/saber/src/saber/bump/${module}.F90#L${i}) | ${purpose} |\n" >> ${autodoc}/${module}.md
               fi
               new_subfunc=false
               type_bound=false
            fi
         fi

         # Reset
         new_purpose=false
      done < ${filename}
   done
done
