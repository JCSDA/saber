#!/bin/bash
#----------------------------------------------------------------------
# Bash shell script: autodoc
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# Directories
src=$1/../../src
doc=$1
mkdir -p ${doc}
autodoc=$1/autodoc
mkdir -p ${autodoc}
dir=(mains/bump saber/bump saber/external saber/gaugrid saber/util)

# Introduction
cat<<EOFINTRO > ${doc}/code_autodoc.md
# Code auto-documentation of Fortran sources

The source directory is organized as follows:
\`\`\`
.
└── saber
    ├── bump:          BUMP core and interfaces
    ├── external:      External tools
    ├── gaugrid:       Gaussian grid tools
    ├── interpolation: Interpolation interface
    ├── oops:          OOPS interface
    └── util:          Shared tools
\`\`\`
and the test directory as follows:
\`\`\`
.
├── mains:     Executables sources
│   └── model: Model-specific interfaces for standalone BUMP
├── testinput: Input YAML files
└── testlist:  List of tests and files used in the CMakeLists.txt
\`\`\`
EOFINTRO

for index in ${!dir[*]}; do
   echo "Directory: "${dir[$index]}
   list=`ls ${src}/${dir[$index]}/*.F90`

   echo -e "## ${dir[$index]}\n" >> ${doc}/code_autodoc.md
   echo -e "| Name | Purpose |" >> ${doc}/code_autodoc.md
   echo -e "| :--: | :---------- |" >> ${doc}/code_autodoc.md

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
            echo "   Module: "${module}
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
               echo -e "| [${module}](autodoc/${module}.md) | ${purpose} |" >> ${doc}/code_autodoc.md

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
                  echo -e "| ${subfunc_type} | [${class}\%] [${subfunc#${class}_}](https://github.com/JCSDA/saber/tree/develop/src/${dir[$index]}/${module}.F90#L${i}) | ${purpose} |" >> ${autodoc}/${module}.md
               else
                  echo -e "| ${subfunc_type} | [${subfunc}](https://github.com/JCSDA/saber/tree/develop/src/${dir[$index]}/${module}.F90#L${i}) | ${purpose} |" >> ${autodoc}/${module}.md
               fi
               new_subfunc=false
               type_bound=false
            fi
         fi

         # Reset
         new_purpose=false
      done < ${filename}
   done
   echo -e "\n" >> ${doc}/code_autodoc.md
done
