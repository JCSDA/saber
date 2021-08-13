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

def sampling_grids(testdata, test, mpi, omp, suffix, testfig):
   """! Plot script for the "sampling grids" files produced by BUMP"""

   # c3 plot
   ic3_plot = [0,2,4,8]

   # c2 plot
   ic2_plot = [0,2,4,8]

   # Loop over files
   first = True
   for impi in range(0, int(mpi)):
      # Open file
      f = Dataset(testdata + "/" + test + "/test_" + mpi + "-" + omp + "_" + suffix + "_" + mpi.zfill(6) + "-" + format(impi+1, '06d') + ".nc", "r", format="NETCDF4")

      # Get _FillValue
      _FillValue = f.__dict__["_FillValue"]

      # Check what is in the file
      if first:
         sc3 = ("lon" in f.variables)
         local = ("lon_local" in f.variables)
         vbal = ("lon_vbal" in f.variables)

      # Get data
      nc0a = f["lon_c0a"].shape[0]
      if sc3:
         nc1a = f["lon"].shape[2]
      if local:
         nc2a = f["lon_local"].shape[1]
      if vbal:
         nc2b = f["lon_vbal"].shape[1]
      if first:
         nl0 = f["gmask_c0a"].shape[0]
         lon_c0 = np.zeros((0))
         lat_c0 = np.zeros((0))
         gmask_c0 = np.zeros((0,nl0))
         if sc3:
            nc3 = f["lon"].shape[1]
            lon = np.zeros((0, nl0, nc3))
            lat = np.zeros((0, nl0, nc3))
         if local:
            nc1max = f["lon_local"].shape[2]
            lon_local = np.zeros((0, nl0, nc1max))
            lat_local = np.zeros((0, nl0, nc1max))
         if vbal:
            nc1max = f["lon_vbal"].shape[2]
            lon_vbal = np.zeros((0, nl0, nc1max))
            lat_vbal = np.zeros((0, nl0, nc1max))
         first = False
      for ic0a in range(0, nc0a):
         lon_c0 = np.append(lon_c0, f["lon_c0a"][ic0a].data)
         lat_c0 = np.append(lat_c0, f["lat_c0a"][ic0a].data)
         gmask_c0 = np.append(gmask_c0, [f["gmask_c0a"][:,ic0a].data], axis=0)
      if sc3:
         for ic1a in range(0, nc1a):
            lon = np.append(lon, [f["lon"][:,:,ic1a].data], axis=0)
            lat = np.append(lat, [f["lat"][:,:,ic1a].data], axis=0)
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
   nc0 = lon_c0.shape[0]
   if sc3:
      nc1 = lon.shape[0]
   if local:
      nc2 = lon_local.shape[0]
   if vbal:
      nc2btot = lon_vbal.shape[0]

   # Set geographical mask
   for ic0 in range(0, nc0):
      gmask_c0[ic0,:] = gmask_c0[ic0,:]*(ic0%2+1)
   gmask_c0 = ma.masked_where(gmask_c0<0.5, gmask_c0, False)

   # Get levels
   levs = f.__dict__["nam_levs"].split(":")

   # Set masked values and levels
   gmask_c0 = ma.masked_invalid(gmask_c0)
   vmax = np.max(gmask_c0)
   gmask_c0 = gmask_c0.filled(fill_value=-1.0e38)
   levels = [-1.0, 0.0, 1.0, 2.0]

   if sc3:
      for ic3 in ic3_plot:
         # Plots
         fig, ax = plt.subplots(nrows=nl0)
         fig.subplots_adjust(hspace=0.4)
         for il0 in range(0, nl0):
            ax[il0].set_title("@ " + str(il0))
            ax[il0].tricontourf(lon_c0, lat_c0, gmask_c0[:,il0], levels=levels, cmap="gray")
            for ic1 in range(0, nc1):
               if (not lon.mask[ic1,il0,ic3]) and (not lat.mask[ic1,il0,ic3]):
                  if (np.abs(lon[ic1,il0,0]-lon[ic1,il0,ic3]) < 180.0):
                     ax[il0].plot([lon[ic1,il0,0], lon[ic1,il0,ic3]], [lat[ic1,il0,0], lat[ic1,il0,ic3]], "k-")
                     ax[il0].plot([lon[ic1,il0,0], lon[ic1,il0,ic3]], [lat[ic1,il0,0], lat[ic1,il0,ic3]], marker="o", color="black", ms=0.2)

         # Save and close figure
         plt.savefig(testfig + "/test_" + mpi + "-" + omp + "_" + suffix + "_c3_" + str(ic3) + ".jpg", format="jpg", dpi=300)
         plt.close()

   if sc3 and (local or vbal):
      # Local sampling

      # Plots
      fig, ax = plt.subplots(nrows=nl0)
      fig.subplots_adjust(hspace=0.4)
      for il0 in range(0, nl0):
         ax[il0].set_title("@ " + str(il0))
         ax[il0].tricontourf(lon_c0, lat_c0, gmask_c0[:,il0], levels=levels, cmap="gray")
         for ic1 in range(0, nc1):
            if (not lon.mask[ic1,il0,0]) and (not lat.mask[ic1,il0,0]):
               ax[il0].plot(lon[ic1,il0,0], lat[ic1,il0,0], marker="o", color="black", ms=0.2)

         if local:
            for ic2 in range(0, nc2):
               if (not lon_local.mask[ic2,il0,0]) and (not lat_local.mask[ic2,il0,0]):
                  ax[il0].plot(lon_local[ic2,il0,0], lat_local[ic2,il0,0], marker="o", color="red", ms=0.2)
         if vbal:
            for ic2 in range(0, nc2btot):
               if (not lon_vbal.mask[ic2,il0,0]) and (not lat_vbal.mask[ic2,il0,0]):
                  ax[il0].plot(lon_vbal[ic2,il0,0], lat_vbal[ic2,il0,0], marker="o", color="red", ms=0.2)

      # Save and close figure
      plt.savefig(testfig + "/test_" + mpi + "-" + omp + "_" + suffix + "_c2_in_c1.jpg", format="jpg", dpi=300)
      plt.close()

   if local:
      # Local mask
      for ic2 in ic2_plot:
         # Plots
         fig, ax = plt.subplots(nrows=nl0)
         fig.subplots_adjust(hspace=0.4)
         for il0 in range(0, nl0):
            ax[il0].set_title("@ " + str(il0))
            ax[il0].tricontourf(lon_c0, lat_c0, gmask_c0[:,il0], levels=levels, cmap="gray")
            for ic1 in range(0, nc1max):
               if (not lon_local.mask[ic2,il0,ic1]) and (not lat_local.mask[ic2,il0,ic1]):
                  if (np.abs(lon_local[ic2,il0,0]-lon_local[ic2,il0,ic1]) < 180.0):
                     ax[il0].plot([lon_local[ic2,il0,0], lon_local[ic2,il0,ic1]], [lat_local[ic2,il0,0], lat_local[ic2,il0,ic1]], "k-")
                     ax[il0].plot([lon_local[ic2,il0,0], lon_local[ic2,il0,ic1]], [lat_local[ic2,il0,0], lat_local[ic2,il0,ic1]], marker="o", color="black", ms=0.2)

         # Save and close figure
         plt.savefig(testfig + "/test_" + mpi + "-" + omp + "_" + suffix + "_c2_" + str(ic2) + ".jpg", format="jpg", dpi=300)
         plt.close()

   if vbal:
      # Vertical balance mask
      for ic2 in ic2_plot:
         # Plots
         fig, ax = plt.subplots(nrows=nl0)
         fig.subplots_adjust(hspace=0.4)
         for il0 in range(0, nl0):
            ax[il0].set_title("@ " + str(il0))
            ax[il0].tricontourf(lon_c0, lat_c0, gmask_c0[:,il0], levels=levels, cmap="gray")
            for ic1 in range(0, nc1max):
               if (not lon_vbal.mask[ic2,il0,ic1]) and (not lat_vbal.mask[ic2,il0,ic1]):
                  if (np.abs(lon_vbal[ic2,il0,0]-lon_vbal[ic2,il0,ic1]) < 180.0):
                     ax[il0].plot([lon_vbal[ic2,il0,0], lon_vbal[ic2,il0,ic1]], [lat_vbal[ic2,il0,0], lat_vbal[ic2,il0,ic1]], "k-")
                     ax[il0].plot([lon_vbal[ic2,il0,0], lon_vbal[ic2,il0,ic1]], [lat_vbal[ic2,il0,0], lat_vbal[ic2,il0,ic1]], marker="o", color="black", ms=0.2)

         # Save and close figure
         plt.savefig(testfig + "/test_" + mpi + "-" + omp + "_" + suffix + "_c2_" + str(ic2) + ".jpg", format="jpg", dpi=300)
         plt.close()
