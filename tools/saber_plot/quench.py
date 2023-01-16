#!/usr/bin/env python3

import argparse
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from netCDF4 import Dataset
import numpy as np
import os
import cartopy.crs as ccrs

def func(args, suffix):
   """! Plot script for the files produced by QUENCH"""

   # Make directory
   if not os.path.exists(args.test):
      os.mkdir(args.test)

   # Open file
   f = Dataset(args.testdata + "/" + args.test + "/" + args.mpi + "-" + args.omp + "_" + suffix + ".nc", "r", format="NETCDF4")

   # Get lon/lat
   lon = f["lon"][:,:]
   lat = f["lat"][:,:]

   # Variables list
   variables = []
   if args.variable is None:
      for var in f.variables:
         if var != "lon" and var != "lat":
            variables.append(var)
   else:
      variables.append(args.variable)

   # Color map
   cmap = matplotlib.cm.get_cmap('coolwarm')

   for var in variables:
      rank = len(f[var].shape)
      if rank == 2:
         # Read variable
         field = f[var][:,:]

         # Min/max values
         vmax = np.max(np.abs(field))
         vmin = -vmax

         # Create figure
         fig,ax = plt.subplots(subplot_kw=dict(projection=ccrs.InterruptedGoodeHomolosine()))
         ax.set_global()
         ax.coastlines(resolution='110m')
         ax.pcolormesh(np.transpose(lon), np.transpose(lat), np.transpose(field), transform=ccrs.PlateCarree(), cmap=cmap, vmin=vmin, vmax=vmax, shading="auto");

         # Set title
         lvmax = np.max(np.abs(field))
         ax.set_title(var + "\n Global amplitude: " + format(vmax, "5.2e") + " \n Level amplitude: " + format(lvmax, "5.2e"))

         # Save and close figure
         plotpath = args.test + "/" + args.mpi + "-" + args.omp + "_" + suffix + "_" + var + ".jpg"
         plt.savefig(plotpath, format="jpg", dpi=300)
         plt.close()
         print(" -> plot produced: " + plotpath)
      else:
         # Read variable
         fields = f[var][:,:,:]

         # Number of levels
         levels = fields.shape[0]

         # Minimum level
         if args.levels_min is None:
            levels_min = 0
         else:
            levels_min = int(args.levels_min)-1

         # Maximum level
         if args.levels_max is None:
            levels_max = levels-1
         else:
            levels_max = int(args.levels_max)-1

         # Step between levels
         if args.levels_step is None:
            levels_step = 1
         else:
            levels_step = int(args.levels_step)

         # Min/max values
         vmax = np.max(np.abs(fields[levels_min:levels_max,:,:]))
         vmin = -vmax

         for jlev in range(max(0,levels_min), min(levels_max, levels), levels_step):
            # Field for this level
            field = fields[jlev,:,:]

            # Create figure
            fig,ax = plt.subplots(subplot_kw=dict(projection=ccrs.InterruptedGoodeHomolosine()))
            ax.set_global()
            ax.coastlines(resolution='110m')
            ax.pcolormesh(np.transpose(lon), np.transpose(lat), np.transpose(field), transform=ccrs.PlateCarree(), cmap=cmap, vmin=vmin, vmax=vmax, shading="auto");

            # Set title
            jlevp1 = str(jlev+1)
            lvmax = np.max(np.abs(field))
            ax.set_title(var + " @ level " + jlevp1 + "\n Global amplitude: " + format(vmax, "5.2e") + " \n Level amplitude: " + format(lvmax, "5.2e"))

            # Save and close figure
            plotpath = args.test + "/" + args.mpi + "-" + args.omp + "_" + suffix + "_" + var + "_" + jlevp1.zfill(3) + ".jpg"
            plt.savefig(plotpath, format="jpg", dpi=300)
            plt.close()
            print(" -> plot produced: " + plotpath)
