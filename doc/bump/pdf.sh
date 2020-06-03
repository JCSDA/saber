#!/bin/bash
#----------------------------------------------------------------------
# Bash script: pdf
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2015-... UCAR, CERFACS, METEO-FRANCE and IRIT
#----------------------------------------------------------------------

if type "pdflatex" > /dev/null ; then
   # Directories
   doc=$1
   latex=$1/latex
   pdf=$1/pdf
   mkdir -p ${pdf}
   cd ${latex}

   for filename in "covariance_filtering" "diffusion_matern_function" "multivariate_localization" "nicas" ; do
      echo "File: "${filename}
      mkdir -p tmp
      pdflatex -output-directory=tmp ${filename}.tex > /dev/null
      bibtex tmp/${filename} > /dev/null
      pdflatex -output-directory=tmp ${filename}.tex > /dev/null
      mv tmp/${filename}.pdf ${pdf}
      rm -fr tmp
   done
else
   echo "pdflatex not found: no pdf documentation"
fi
