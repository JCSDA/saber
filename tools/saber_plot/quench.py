#!/usr/bin/env python3

import argparse
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import os
import cartopy.crs as ccrs

def func(args):
   """! Plot script for the files produced by QUENCH"""

   # Open longitudes
   with open(args.filepath + "_lon") as f:
      data = f.read()
      nodes = int(data.split(",", 1)[0].replace("size: ", ""))
      lonstring = data.split(",", 1)[1].replace("values: [ ", "").replace(" ]", "")
      lon = np.zeros((nodes))
      inode = 0
      for item in lonstring.split():
         lon[inode] = float(item)
         inode += 1

   # Open latitudes
   with open(args.filepath + "_lat") as f:
      data = f.read()
      nodes = int(data.split(",", 1)[0].replace("size: ", ""))
      latstring = data.split(",", 1)[1].replace("values: [ ", "").replace(" ]", "")
      lat = np.zeros((nodes))
      inode = 0
      for item in latstring.split():
         lat[inode] = float(item)
         inode += 1

   # Open variable
   with open(args.filepath + "_" + args.variable) as f:
      data = f.read()
      n = int(data.split(",", 1)[0].replace("size: ", ""))
      levels = int(n/nodes)
      varstring = data.split(",", 1)[1].replace("values: [ ", "").replace(" ]", "")
      field = np.zeros((nodes,levels))
      inode = 0
      ilevel = 0
      for item in varstring.split():
         field[inode,ilevel] = float(item)
         ilevel += 1
         if ilevel == levels:
            ilevel = 0
            inode +=1

   # Color map
   cmap = matplotlib.cm.get_cmap('coolwarm')

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
      ax.set_title(args.variable + " @ level " + jlevp1 + "\n Global amplitude: " + format(vmax, "5.2e") + " \n Level amplitude: " + format(lvmax, "5.2e"))

      # Save and close figure
      if args.output is None:
         plotpath = os.path.splitext(os.path.basename(args.filepath))[0] + "_" + args.variable + "_" + jlevp1.zfill(3) + ".jpg"
      else:
         plotpath = args.output + "_" + args.variable + "_" + jlevp1.zfill(3) + ".jpg"
      plt.savefig(str(plotpath), format="jpg", dpi=300)
      plt.close()
      print(" -> plot produced: " + plotpath)
