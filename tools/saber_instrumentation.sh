#!/usr/bin/env bash
#----------------------------------------------------------------------
# Bash script: saber_instrumentation
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

# To run the script over all the source files:
# for file in `find .. -name '*.F90' -o -name '*.fypp'` ; do echo $file ; ./saber_instrumentation.sh $file ; done
#
# To get the list of all subroutines/functions in the pre-processed files:
# grep -si probe%in\( /mnt/data/build/gnu_9.3.0/bundle/saber/src/saber/*/*.F90 | gawk 'BEGIN{FS="'\''"}{print $2}' > subr_list.txt ; sed -i 's/^/#:set subr_list = subr_list + [\"/' subr_list.txt; sed -i 's/$/\"]/' subr_list.txt


# File to instrument
filename=$1

# Nothing to do for these files
if [[ ${filename} == *"/tools_kinds.fypp" ]] || [[ ${filename} == *"/type_probe.fypp" ]] || [[ ${filename} == *"/instrumentation"*"fypp" ]] || [[ ${filename} == *"/generics.fypp" ]] ; then
   exit
fi

# First pass
mv ${filename} ${filename}.old
touch ${filename}
OLDIFS=$IFS
IFS=''
cat ${filename}.old | while read -r line ; do
   if [[ ${line} == "subroutine "* ]] || [[ ${line} == "function "* ]] ; then
      tmp=${line#* }
      subr=${tmp%%(*}
   fi

   # Remove subr argument for mpl_nc_* subroutines/functions
   if [[ ${line} == *"mpl_nc_"* ]] || [[ ${line} == *"mpl%nc_"* ]] ; then
      line=${line/"subr,"/""}
   fi

   if [[ ${subr} != "mpl_abort" ]] && [[ ${subr} != "mpl_warning" ]] && [[ ${subr} != "mpl_ncerr" ]] ; then
      # Replace the subr argument
      if [[ ${line} == *"subr,"* ]] ; then
         line=${line/"subr,"/"'\${subr}\$',"}
         if [[ ${#line} -gt 132 ]] ; then
            echo ${line%%"$"*}"\${subr}\$', &" >> ${filename}
            line=" & "${line##*$\',}
         fi
      fi
   fi

   # Remove trailing whitespace
   line="${line%"${line##*[![:space:]]}"}"

   # Remove subr variable
   if [[ ${subr} == "mpl_abort" ]] || [[ ${subr} == "mpl_warning" ]] || [[ ${subr} == "mpl_ncerr" ]] ; then
      # Print line
      echo ${line} >> ${filename}
   else
      if [[ ${line} != *"character"*"subr"* ]] ; then
         # Print line
         echo ${line} >> ${filename}
      fi
   fi
done
IFS=$OLDIFS
rm -f ${filename}.old

# Second pass
mv ${filename} ${filename}.old
touch ${filename}
OLDIFS=$IFS
IFS=''
status=none
cat ${filename}.old | while read -r line ; do
   # Remove "! Local variable" when only subr was declared
   if test "${status}" = "wait" ; then
      status=none
   fi
   if test "${status}" = "local_variables" ; then
      if [[ ${line} == "" ]] ; then
         status=wait
      else
         echo ${line_save} >> ${filename}
         status=none
      fi
   fi
   if [[ ${line} == "! Local variable"* ]] ; then
      status="local_variables"
      line_save=${line}
   fi

   if test "${status}" = "none" ; then
      # Print line
      echo ${line} >> ${filename}
   fi
done
IFS=$OLDIFS
rm -f ${filename}.old

# Nothing to do for these files
if [[ ${filename} == *".F90" ]] ; then
   exit
fi

# Third pass
mv ${filename} ${filename}.old
touch ${filename}
OLDIFS=$IFS
IFS=''
new_subr=false
add_module=false
module_old="aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
start_modules=false
status=none
line_prev=""
init=true

cat ${filename}.old | while read -r line ; do
   # Add first line
   if test "${init}" = "true" ; then
      if [[ ${line} == "!----------------------------------------------------------------------" ]] ; then
         echo "#:include '../instrumentation.fypp'" >> ${filename}
         init=false
      else
         if [[ ${line} == "#:include '../instrumentation.fypp'" ]] ; then
            init=false
         fi
      fi
   fi

   # Add type_probe module
   if [[ ${line:0:3} == "use" ]] ; then
      start_modules=true
      module=`echo ${line} | gawk '{print $2}'`
      if [[ ${module} < ${module_old} ]] ; then
         echo "Error: modules are not ordered alphabetically => ${module} < ${module_old}"
         exit 1
      fi
      module_old=${module}
      if [[ ${module} == "type_probe," ]] ; then
         add_module=true
      fi
      if [[ ${add_module} == "false" ]] ; then
         if [[ ${module} > "type_probe," ]] ; then
            echo "use type_probe, only: probe" >> ${filename}
            add_module=true
            module_old="type_probe,"
         fi
      fi
   else
      if [[ ${start_modules}  == "true" ]] && [[ ${add_module} == "false" ]] && [[ ${line} != " & "* ]] && [[ ${line} != "!$ "* ]] ; then
         echo "use type_probe, only: probe" >> ${filename}
         add_module=true
      fi
   fi

   if [[ ${line} == "subroutine "* ]] || [[ ${line} == "function "* ]] ; then
      tmp=${line#* }
      subr=${tmp%%(*}
      new_subr=true
      new_declaration=false
      new_declaration_prev=false
      instrumentation_needed=false
      status=none
   fi

   if [[ ${subr} != "nam_read" ]] && [[ ${subr} != "mpl_timings" ]] && [[ ${subr} != "mpl_warning" ]] && [[ ${subr} != "mpl_abort" ]] && [[ ${subr} != "mpl_ncerr" ]] ; then
      if test "${status}" = "returned_variable_2"; then
         # Nothing to do (returned variable)
         status="possible_end_of_declarations"
      fi

      if test "${status}" = "returned_variable_1"; then
         # Nothing to do (returned variable)
         if [[ ${line} != *"&"* ]] && [[ ${line} != "#"* ]]  ; then
            status="returned_variable_2"
         fi
      fi

      if test "${status}" = "possible_end_of_declarations"; then
         if [[ ${line} == "! Returned variable" ]] ; then
            status="returned_variable_1"
         else
            if [[ ${line} == "! Local variable"* ]] ; then
               status="local_variables"
            else
               # Instrumentation
               status="instrumentation_in"
            fi
         fi
      fi

      if test "${status}" != "returned_variable_1"; then
         if test "${new_subr}" = "true"; then
            if [[ ${line} == "integer"* ]] || [[ ${line} == "real"* ]] || [[ ${line} == "logical"* ]] || [[ ${line} == "character"* ]] || [[ ${line} == "type"* ]] || [[ ${line} == "class"* ]] || [[ ${line} == *"ftype"* ]] || [[ ${line} == "#"*"integer"* ]] || [[ ${line} == "#"*"real"* ]] || [[ ${line} == "#"*"logical"* ]] || [[ ${line} == "#"*"character"* ]] || [[ ${line} == "#"*"type"* ]]; then
               new_declaration=true
            else
               new_declaration=false
            fi
            if test "${new_declaration}" = "false" && "${new_declaration_prev}" = "true" ; then
               status="possible_end_of_declarations"
            fi
         fi

         if test "${new_declaration}" = "true"; then
            new_declaration_prev=true
         else
            new_declaration_prev=false
         fi
      fi

      if test "${status}" = "local_variables"; then
         if [[ ${line} == "integer"* ]] || [[ ${line} == "real"* ]] || [[ ${line} == "logical"* ]] || [[ ${line} == "character"* ]] || [[ ${line} == "type"* ]] || [[ ${line} == "class"* ]] || [[ ${line} == *"ftype"* ]] || [[ ${line} == "#"*"integer"* ]] || [[ ${line} == "#"*"real"* ]] || [[ ${line} == "#"*"logical"* ]] || [[ ${line} == "#"*"character"* ]] || [[ ${line} == "#"*"type"* ]]; then
            # Nothing to do
            status="local_variables"
         else
            # Instrumentation
            status="instrumentation_in_next"
         fi
      fi

      if test "${status}" = "instrumentation_in"; then
         if [[ ${line} != "! Set name" ]] ; then
            echo "! Set name" >> ${filename}
            echo "#:set subr = '${subr}'" >> ${filename}
            echo "" >> ${filename}
            if [[ ${subr} == "bint_"* ]] ; then
               echo "! Get instance index" >> ${filename}
               echo "call probe%get_instance(bint%bump%iinst)" >> ${filename}
               echo "" >> ${filename}
            fi
            if [[ ${subr} == "bump_"* ]] ; then
               echo "! Get instance index" >> ${filename}
               echo "call probe%get_instance(bump%iinst)" >> ${filename}
               echo "" >> ${filename}
            fi
            if [[ ${subr} == "gaugrid_"* ]] ; then
               echo "! Get instance index" >> ${filename}
               echo "call probe%get_instance(self%iinst)" >> ${filename}
               echo "" >> ${filename}
            fi
            echo "! Probe in" >> ${filename}
            echo "@:probe_in()" >> ${filename}
            echo "" >> ${filename}
            instrumentation_needed=true
         fi
         status="done"
         new_subr=0
         new_declaration=false
         new_declaration_prev=false
      fi

      if [[ ${line} == "end subroutine "* ]] || [[ ${line} == "end function "* ]] ; then
         if test "${instrumentation_needed}" = "true" ; then
            # Instrumentation
            echo "! Probe out" >> ${filename}
            echo "@:probe_out()" >> ${filename}
            echo "" >> ${filename}
         fi
      fi
   fi

   if [[ ${line} == *"return" ]] ; then
      if [[ "${instrumentation_needed}" == "true" ]] && [[ ${line_prev} != "@:probe_out()" ]] ; then
         # Instrumentation
         echo "@:probe_out()" >> ${filename}
      fi
   fi

   # Print line
   echo ${line} >> ${filename}

   # Save line
   line_prev=${line}
done
IFS=$OLDIFS
rm -f ${filename}.old
