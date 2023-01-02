#!/usr/bin/env python3

import argparse
from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import numpy.ma as ma
import os
import cartopy.crs as ccrs

def bump_sampling_grids(args, suffix):
   """! Plot script for the "sampling grids" files produced by BUMP"""

   # Make directory
   if not os.path.exists(args.test):
      os.mkdir(args.test)

   # Define colors for angular sectors
   colors = []
   for item in mcolors.TABLEAU_COLORS:
      colors.append(item)

   # il0
   il0 = 0

   # nc1 plot
   nc1_plot = 10

   # nc2 plot
   nc2_plot = 10

   # Loop over files
   first = True
   for impi in range(0, int(args.mpi)):
      # Open file
      f = Dataset(args.testdata + "/" + args.test + "/" + args.mpi + "-" + args.omp + "_" + suffix + "_local_" + args.mpi.zfill(6) + "-" + format(impi+1, '06d') + ".nc", "r", format="NETCDF4")

      # Get _FillValue
      _FillValue = f.__dict__["_FillValue"]

      # Check what is in the file
      if first:
         sc3 = ("lon" in f.variables)
         local = ("lon_local" in f.variables)
         vbal = ("lon_vbal" in f.variables)

      # Get data
      if sc3:
         nc1a = f["lon"].shape[3]
      if local:
         nc2a = f["lon_local"].shape[1]
      if vbal:
         nc2b = f["lon_vbal"].shape[1]
      if first:
         nl0 = f["lon"].shape[0]
         if sc3:
            nc3 = f["lon"].shape[2]
            nc4 = f["lon"].shape[1]
            lon = np.zeros((0, nl0, nc4, nc3))
            lat = np.zeros((0, nl0, nc4, nc3))
         if local:
            nc1max = f["lon_local"].shape[2]
            lon_local = np.zeros((0, nl0, nc1max))
            lat_local = np.zeros((0, nl0, nc1max))
         if vbal:
            nc1max = f["lon_vbal"].shape[2]
            lon_vbal = np.zeros((0, nl0, nc1max))
            lat_vbal = np.zeros((0, nl0, nc1max))
         first = False
      if sc3:
         for ic1a in range(0, nc1a):
            lon = np.append(lon, [f["lon"][:,:,:,ic1a].data], axis=0)
            lat = np.append(lat, [f["lat"][:,:,:,ic1a].data], axis=0)
      if local:
         for ic2a in range(0, nc2a):
            lon_local = np.append(lon_local, [f["lon_local"][:,ic2a,:].data], axis=0)
            lat_local = np.append(lat_local, [f["lat_local"][:,ic2a,:].data], axis=0)
      if vbal:
         for ic2b in range(0, nc2b):
            lon_vbal = np.append(lon_vbal, [f["lon_vbal"][:,ic2b,:].data], axis=0)
            lat_vbal = np.append(lat_vbal, [f["lat_vbal"][:,ic2b,:].data], axis=0)

   # Apply mask
   if sc3:
      lon = ma.masked_where(lon==_FillValue, lon, False)
      lat = ma.masked_where(lat==_FillValue, lat, False)
   if local:
      lon_local = ma.masked_where(lon_local==_FillValue, lon_local, False)
      lat_local = ma.masked_where(lat_local==_FillValue, lat_local, False)
   if vbal:
      lon_vbal = ma.masked_where(lon_vbal==_FillValue, lon_vbal, False)
      lat_vbal = ma.masked_where(lat_vbal==_FillValue, lat_vbal, False)

   # Total sizes
   if sc3:
      nc1 = lon.shape[0]
   if local:
      nc2 = lon_local.shape[0]
   if vbal:
      nc2btot = lon_vbal.shape[0]

   # To plot
   if sc3:
      ic1_plot = np.linspace(0, nc1-1, num=nc1_plot).astype(int)
   if local:
      ic2_plot = np.linspace(0, nc2-1, num=nc2_plot).astype(int)
   if vbal:
      ic2btot_plot = np.linspace(0, nc2btot-1, num=nc2_plot).astype(int)

   if sc3:
      for ic1 in ic1_plot:
         # Plots
         fig, ax = plt.subplots(subplot_kw=dict(projection=ccrs.NearsidePerspective(central_longitude=lon[ic1,0,0,0], central_latitude=lat[ic1,0,0,0])))
         ax.set_title("Level " + str(il0+1))
         ax.set_global()
         ax.coastlines(resolution='110m')
         for ic4 in range(nc4):
            for ic3 in range(nc3):
               if (not lon.mask[ic1,il0,ic4,ic3]) and (not lat.mask[ic1,il0,ic4,ic3]):
                  ax.plot([lon[ic1,il0,0,0], lon[ic1,il0,ic4,ic3]], [lat[ic1,il0,0,0], lat[ic1,il0,ic4,ic3]], colors[ic4], transform=ccrs.Geodetic())

         # Save and close figure
         plotpath = args.test + "/" + args.mpi + "-" + args.omp + "_" + suffix + "_c1_" + str(ic1) + ".jpg"
         plt.savefig(plotpath, format="jpg", dpi=300)
         plt.close()
         print(" -> plot produced: " + plotpath)

   if sc3 and (local or vbal):
      # Local sampling

      # Plots
      fig, ax = plt.subplots(subplot_kw=dict(projection=ccrs.NearsidePerspective()))
      ax.set_title("Level " + str(il0+1))
      ax.set_global()
      ax.coastlines(resolution='110m')
      for ic1 in range(0, nc1):
         if (not lon.mask[ic1,il0,0,0]) and (not lat.mask[ic1,il0,0,0]):
            ax.plot(lon[ic1,il0,0,0], lat[ic1,il0,0,0], marker="o", color="black", ms=0.2, transform=ccrs.Geodetic())

      if local:
         for ic2 in range(0, nc2):
            if (not lon_local.mask[ic2,il0,0]) and (not lat_local.mask[ic2,il0,0]):
               ax.plot(lon_local[ic2,il0,0], lat_local[ic2,il0,0], marker="o", color="red", ms=0.2, transform=ccrs.Geodetic())
      if vbal:
         for ic2 in range(0, nc2btot):
            if (not lon_vbal.mask[ic2,il0,0]) and (not lat_vbal.mask[ic2,il0,0]):
               ax.plot(lon_vbal[ic2,il0,0], lat_vbal[ic2,il0,0], marker="o", color="red", ms=0.4, transform=ccrs.Geodetic())

      # Save and close figure
      plotpath = args.test + "/" + args.mpi + "-" + args.omp + "_" + suffix + "_c2_in_c1.jpg"
      plt.savefig(plotpath, format="jpg", dpi=300)
      plt.close()
      print(" -> plot produced: " + plotpath)

   if local:
      # Local mask
      for ic2 in ic2_plot:
         # Plots
         fig, ax = plt.subplots(subplot_kw=dict(projection=ccrs.NearsidePerspective(central_longitude=lon_local[ic2,il0,0], central_latitude=lat_local[ic2,il0,0])))
         ax.set_title("Level " + str(il0+1))
         ax.set_global()
         ax.coastlines(resolution='110m')
         for ic1 in range(0, nc1max):
             if (not lon_local.mask[ic2,il0,ic1]) and (not lat_local.mask[ic2,il0,ic1]):
               ax.plot(lon_local[ic2,il0,ic1], lat_local[ic2,il0,ic1], marker="o", color="black", ms=0.2, transform=ccrs.Geodetic())
         ax.plot(lon_local[ic2,il0,0], lat_local[ic2,il0,0], marker="o", color="red", ms=0.4, transform=ccrs.Geodetic())

         # Save and close figure
         plotpath = args.test + "/" + args.mpi + "-" + args.omp + "_" + suffix + "_c2_" + str(ic2) + ".jpg"
         plt.savefig(plotpath, format="jpg", dpi=300)
         plt.close()
         print(" -> plot produced: " + plotpath)

   if vbal:
      # Vertical balance mask
      for ic2 in ic2btot_plot:
         # Plots
         fig, ax = plt.subplots(subplot_kw=dict(projection=ccrs.NearsidePerspective(central_longitude=lon_vbal[ic2,il0,0], central_latitude=lat_vbal[ic2,il0,0])))
         ax.set_title("Level " + str(il0+1))
         ax.set_global()
         ax.coastlines(resolution='110m')
         for ic1 in range(0, nc1max):
            if (not lon_vbal.mask[ic2,il0,ic1]) and (not lat_vbal.mask[ic2,il0,ic1]):
               ax.plot(lon_vbal[ic2,il0,ic1], lat_vbal[ic2,il0,ic1], marker="o", color="black", ms=0.2, transform=ccrs.Geodetic())
         ax.plot(lon_vbal[ic2,il0,0], lat_vbal[ic2,il0,0], marker="o", color="red", ms=0.4, transform=ccrs.Geodetic())

         # Save and close figure
         plotpath = args.test + "/" + args.mpi + "-" + args.omp + "_" + suffix + "_c2_" + str(ic2) + ".jpg"
         plt.savefig(plotpath, format="jpg", dpi=300)
         plt.close()
         print(" -> plot produced: " + plotpath)
