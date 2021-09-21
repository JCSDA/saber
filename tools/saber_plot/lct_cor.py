#!/usr/bin/env python3

import argparse
from netCDF4 import Dataset
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import numpy.ma as ma
import os

def lct_cor(testdata, test, mpi, omp, suffix, testfig):
   """! Plot script for the "LCT correlation" files produced by BUMP"""

   # Processor, index and level to plot
   iproc_plot = 1
   ic1a_plot = 1
   il0_plot = 1

   # Open file
   f = Dataset(testdata + "/" + test + "/test_" + mpi + "-" + omp + "_" + suffix + "_" + mpi.zfill(6) + "-" + format(iproc_plot, '06d') + ".nc", "r", format="NETCDF4")

   # Get _FillValue
   _FillValue = f.__dict__["_FillValue"]

   for group in f.groups:
      # Get lon/lat/levels
      lon = f.groups[group]["lon"][:,:]
      lat = f.groups[group]["lat"][:,:]
      l0rl0_to_l0 = f.groups[group]["l0rl0_to_l0"][:,:]

      # Get vertical unit
      vunit = f.groups[group]["vunit"][:]

      # Get number of levels
      nl0 = vunit.shape[0]

      for fieldname in ["raw", "fit", "fit_filt"]:
         if fieldname in f.groups[group].variables:
            # Get field
            field = f.groups[group][fieldname][:,:,:,:]

            # Get dimensions
            nl0r = field.shape[2]
            nc3 = field.shape[3]

            # Set masked values and levels
            field = ma.masked_invalid(field)
            vmax = np.max(np.abs(field))
            if (vmax > 0.0):
               levels = np.linspace(-vmax, vmax, 21)
            else:
               levels = np.linspace(-1.0, 1.0, 3)
            field = field.filled(fill_value=-1.0e38)

            # Plots
            fig, ax = plt.subplots(nrows=nl0r)
            fig.subplots_adjust(hspace=0.4, right=0.8)
            for il0r in range(0, nl0r):
               ax[il0r].set_title(group + " - " + fieldname + " @ " + str(l0rl0_to_l0[il0_plot,il0r]))
               im = ax[il0r].tricontourf(lon[ic1a_plot,:], lat[ic1a_plot,:], field[il0_plot,ic1a_plot,il0r,:], levels=levels, cmap="bwr")

            # Colorbar
            cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
            fig.colorbar(im, cax=cbar_ax)

            # Save and close figure
            plt.savefig(testfig + "/test_" + mpi + "-" + omp + "_" + suffix + "_" + group + "_" + fieldname + ".jpg", format="jpg", dpi=300)
            plt.close()
