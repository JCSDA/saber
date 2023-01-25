#!/bin/bash

# Loop over files
for file in `ls *.yaml *.yml`; do
#  file=error_covariance_training_bump_nicas_5.yaml
  echo "File: "$file
  # Rename bump blocks in alias items
  for pattern in `grep -si "in code:" $file | sed 's/.*://'`; do
    new_pattern=`echo $pattern | sed 's/-.*//'`
    sed -i -e s/'${pattern}'/'${new_pattern}'/g $file
  done

  # Rename multivariate strategy
  sed -i -e s/'multivariate strategy: specific_univariate'/'multivariate strategy: univariate'/g $file
  sed -i -e s/'multivariate strategy: specific_multivariate'/'multivariate strategy: multivariate'/g $file
  sed -i -e s/'multivariate strategy: common_weighted'/'multivariate strategy: duplicated and weighted'/g $file
  sed -i -e s/'multivariate strategy: common'/'multivariate strategy: duplicated'/g $file

  # Rename variables into groups
  foundit=false
  mv $file $file.bak
  IFS=''
  while read line; do
    if test "${foundit}" = "true"; then
      new_indentation=`echo $line | awk -F"[ ]" '{for(i=1;i<=NF && ($i=="");i++);print i-1}'`
      if [[ ${new_indentation} -eq ${indentation} ]]; then
        # Same indentation
        if [[ ${line:$indentation:1} == "-" ]]; then
          echo $line | sed -e s/"- variables:"/"- groups:"/g >> $file
        else
          echo $line >> $file
          foundit=false
        fi
      elif [[ ${new_indentation} -gt ${indentation} ]]; then
        # Larger indentation
        echo $line >> $file
      elif [[ ${new_indentation} -lt ${indentation} ]]; then
        # Lower indentation
        echo $line >> $file
        foundit=false
      fi
    else
      echo $line >> $file
    fi
    if [[ $line =~ "horizontal length-scale:" ]]; then
      foundit=true
      indentation=`echo $line | awk -F"[ ]" '{for(i=1;i<=NF && ($i=="");i++);print i-1}'`
    fi
    if [[ $line =~ "vertical length-scale:" ]]; then
      foundit=true
      indentation=`echo $line | awk -F"[ ]" '{for(i=1;i<=NF && ($i=="");i++);print i-1}'`
    fi
    if [[ $line =~ "interpolation type:" ]]; then
      foundit=true
      indentation=`echo $line | awk -F"[ ]" '{for(i=1;i<=NF && ($i=="");i++);print i-1}'`
    fi
    if [[ $line =~ "minimum level:" ]]; then
      foundit=true
      indentation=`echo $line | awk -F"[ ]" '{for(i=1;i<=NF && ($i=="");i++);print i-1}'`
    fi
    if [[ $line =~ "maximum level:" ]]; then
      foundit=true
      indentation=`echo $line | awk -F"[ ]" '{for(i=1;i<=NF && ($i=="");i++);print i-1}'`
    fi
  done < $file.bak
  if test "$(tail -c 1 $file.bak)"; then
    tail -n 1 $file.bak >> $file
  fi
  rm -f $file.bak
done