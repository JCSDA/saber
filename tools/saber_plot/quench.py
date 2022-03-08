#!/usr/bin/env python3

import argparse
from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import os
import cartopy.crs as ccrs

def func(args):
   """! Plot script for the files produced by QUENCH"""

   # Open file
   f = Dataset(args.filepath, "r", format="NETCDF4")

   # Get lon/lat
   lon = f["lon"][:]
   lat = f["lat"][:]

   # Color map
   cmap = matplotlib.cm.get_cmap('coolwarm')

   # Minimum level
   if args.levels_min is None:
      levels_min = 0
   else:
      levels_min = int(args.levels_min)-1

   # Maximum level
   if args.levels_max is None:
      levels_max = 1000
   else:
      levels_max = int(args.levels_max)-1

   # Step between levels
   if args.levels_step is None:
      levels_step = 1
   else:
      levels_step = int(args.levels_step)

   for var in f.variables:
      if (var != "lon" and var != "lat"):
         # Read variable
         field = f[var][:,:]

         # Number of nodes and levels
         nodes = field.shape[0]
         levels = field.shape[1]

         # Min/max values
         vmax = np.max(np.abs(field[:, levels_min:levels_max]))
         vmin = -vmax

         for jlev in range(max(0,levels_min), min(levels_max, levels), levels_step):
            # Create figure
            fig,ax = plt.subplots(subplot_kw=dict(projection=ccrs.InterruptedGoodeHomolosine()))
            ax.set_global()
            ax.coastlines(resolution='110m')
            for jnode in range(0, nodes):
               col = cmap((field[jnode, jlev]-vmin)/(vmax-vmin))
               ax.plot(lon[jnode], lat[jnode], color=col, marker='.', markersize=1.0, transform=ccrs.PlateCarree())

            # Set title
            jlevp1 = str(jlev+1)
            lvmax = np.max(np.abs(field[:, jlev]))
            ax.set_title(var + " @ level " + jlevp1 + "\n Global amplitude: " + format(vmax, "5.2e") + " \n Level amplitude: " + format(lvmax, "5.2e"))

            # Save and close figure
            if args.output is None:
               plotpath = os.path.splitext(os.path.basename(args.filepath))[0] + "_" + var + "_" + jlevp1.zfill(3) + ".jpg"
            else:
               plotpath = args.output + "_" + var + "_" + jlevp1.zfill(3) + ".jpg"
            plt.savefig(str(plotpath), format="jpg", dpi=300)
            plt.close()
            print(" -> plot produced: " + plotpath)
