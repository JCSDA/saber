#!/usr/bin/env python3

import argparse
from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import os

def diag(testdata, test, mpi, omp, suffix, testfig):
   """! Plot script for the "diagnostic" files produced by BUMP"""

   # Open file
   f = Dataset(testdata + "/" + test + "/test_" + mpi + "-" + omp + "_" + suffix + ".nc", "r", format="NETCDF4")

   # Get _FillValue
   _FillValue = f.__dict__["_FillValue"]

   # Get vertical unit
   vunit = f["vunit"][:]
   vunitmin = np.min(vunit)
   vunitmax = np.max(vunit)

   # Get number of levels
   nl0 = vunit.shape[0]

   # Profiles only
   if nl0 == 1:
      return

   # Diagnostics list
   diag_list = ["coef_ens","fit_rh","fit_rv"]   

   # Plots
   for diag in diag_list:
      fig, ax = plt.subplots()
      fig.subplots_adjust(right=0.8)
      cmap = matplotlib.cm.get_cmap('Spectral')
      ax.set_title(diag)
      ax.set_ylim([vunitmin,vunitmax])
      valid = False
      for group in f.groups:
         for subgroup in f.groups[group].groups:
            if (diag in f.groups[group].groups[subgroup].variables):
               ax.plot(f.groups[group].groups[subgroup][diag][:], vunit, label=group + " - " + subgroup)
               valid = True

      if (valid):
         # Single legend
         handles, labels = ax.get_legend_handles_labels()
         fig.legend(handles, labels, loc='center right')

         # Save and close figure
         plt.savefig(testfig + "/test_" + mpi + "-" + omp + "_" + suffix + "_" + diag + ".jpg", format="jpg", dpi=300)
         plt.close()
      else:
         # Just close the figure
         plt.close()
