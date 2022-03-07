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

def bump_contour_centered(args, suffix):
   """! Plot script for the "centered countour" files produced by BUMP"""

   # Open file
   f = Dataset(args.testdata + "/" + args.test + "/test_" + args.mpi + "-" + args.omp + "_" + suffix + ".nc", "r", format="NETCDF4")

   # Get lon/lat
   lon = f["lon"][:]
   lat = f["lat"][:]

   # Get vertical unit
   vunit = f["vunit"][:,:]

   # Get number of levels
   nl0 = vunit.shape[0]

   for group in f.groups:
      for var in f.groups[group].variables:
         # Read variable
         field = f.groups[group][var][:,:]

         # Set masked values and levels
         field = ma.masked_invalid(field)
         vmax = np.max(np.abs(field))
         if (vmax > 0.0):
            levels = np.linspace(-vmax, vmax, 21)
         else:
            levels = np.linspace(-1.0, 1.0, 3)
         field = field.filled(fill_value=-1.0e38)

         # Plots
         fig, ax = plt.subplots(nrows=nl0)
         fig.subplots_adjust(hspace=0.4, right=0.8)
         for il0 in range(0, nl0):
            ax[il0].set_title(group + " - " + var + " @ " + str(il0))
            im = ax[il0].tricontourf(lon, lat, field[il0,:], levels=levels, cmap="bwr")

         # Colorbar
         cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
         fig.colorbar(im, cax=cbar_ax)

         # Save and close figure
         plotpath = args.test + "_" + args.mpi + "-" + args.omp + "_" + suffix + "_" + group + "_" + var + ".jpg"
         plt.savefig(plotpath, format="jpg", dpi=300)
         plt.close()

      for subgroup in f.groups[group].groups:
         for var in f.groups[group].groups[subgroup].variables:
            # Read variable
            field = f.groups[group].groups[subgroup][var][:,:]

            # Set masked values and levels
            field = ma.masked_invalid(field)
            vmax = np.max(np.abs(field))

            levels = np.linspace(-vmax, vmax, 21)
            field = field.filled(fill_value=-1.0e38)

            # Plots
            fig, ax = plt.subplots(nrows=nl0)
            fig.subplots_adjust(hspace=0.4, right=0.8)
            for il0 in range(0, nl0):
               ax[il0].set_title(group + " - " + subgroup + " - " + var + " @ " + str(il0))
               im = ax[il0].tricontourf(lon, lat, field[il0,:], levels=levels, cmap="bwr")

            # Colorbar
            cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
            fig.colorbar(im, cax=cbar_ax)

            # Save and close figure
            plotpath = args.test + "_" + args.mpi + "-" + args.omp + "_" + suffix + "_" + group + "_" + subgroup + "_" + var + ".jpg"
            plt.savefig(plotpath, format="jpg", dpi=300)
            plt.close()
