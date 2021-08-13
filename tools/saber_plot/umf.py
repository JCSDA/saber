#!/usr/bin/env python3

import argparse
from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import os

def umf(testdata, test, mpi, omp, suffix, testfig):
   """! Plot script for the "univariate moments fields"  files produced by BUMP"""

   # Open file
   f = Dataset(testdata + "/" + test + "/test_" + mpi + "-" + omp + "_" + suffix + ".nc", "r", format="NETCDF4")

   # Get _FillValue
   _FillValue = f.__dict__["_FillValue"]

   # Get lon/lat
   lon = f["lon"][:]
   lat = f["lat"][:]

   # Variables
   var_list = ["m2","m4","kurt"]

   for group in f.groups:
      # Get number of levels
      nl0 = f.groups[group]["vunit"][:,:].shape[0]

      for var in f.groups[group].variables:
         # Read variable
         field = f.groups[group][var][:,:]

         # Set masked values and levels
         field = ma.masked_invalid(field)
         vmax = np.max(field)
         levels = np.linspace(0, vmax, 21)
         field = field.filled(fill_value=-1.0e38)

         # Plots
         fig, ax = plt.subplots(nrows=nl0)
         fig.subplots_adjust(hspace=0.4, right=0.8)
         for il0 in range(0, nl0):
            ax[il0].set_title(var + " @ " + str(il0))
            im = ax[il0].tricontourf(lon, lat, field[il0,:], levels=levels, cmap="YlOrRd")

         # Colorbar
         cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
         fig.colorbar(im, cax=cbar_ax)

         # Save and close figure
         plt.savefig(testfig + "/test_" + mpi + "-" + omp + "_" + suffix + "_" + var + ".jpg", format="jpg", dpi=300)
         plt.close()
