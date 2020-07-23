#!/usr/bin/env python3

import argparse
from netCDF4 import Dataset
import Ngl
import numpy as np
import os

def randomization(testdata, test, mpi, omp, suffix, testfig):
   # Open file
   f = Dataset(testdata + "/" + test + "/test_" + mpi + "-" + omp + "_" + suffix + ".nc", "r", format="NETCDF4")

   # Get factors
   nefac = f["nefac"][:]

   # Get number of tests
   nfac = nefac.shape[0]

   # Get number of factors
   nfac = nefac.shape[0]
   x = Ngl.fspan(1, nfac, nfac)

   # Get data, compute mean and standard deviation
   mse_mean = np.zeros((0, nfac))
   mse_mean = np.append(mse_mean, [np.mean(f["mse"][:,:], axis=1)], axis=0)
   mse_mean = np.append(mse_mean, [np.mean(f["mse_th"][:,:], axis=1)], axis=0)
   mse_std = np.zeros((0, nfac))
   mse_std = np.append(mse_std, [np.std(f["mse"][:,:], axis=1)], axis=0)
   mse_std = np.append(mse_std, [np.std(f["mse_th"][:,:], axis=1)], axis=0)

   # XY resources
   xyres = Ngl.Resources()
   xyres.nglFrame = False
   xyres.nglDraw = False
   xyres.xyMarkLineMode = "Markers"
   xyres.xyMarker = 1
   xyres.xyMarkerColors = ["black", "red"]
   xyres.xyMarkerSizeF = 0.05
   xyres.tiXAxisString = "Number of members"
   xyres.tiYAxisString = "Mean Square Error"
   xyres.nglYRefLine = 1.0
   xyres.trXMinF = 0.5
   xyres.trXMaxF = float(nfac)+0.5
   xyres.tmXBMode = "Explicit"
   xyres.tmXBValues = x
   xyres.tmXBLabels = nefac

   # Polyline resources
   plres = Ngl.Resources()
   plres.gsLineThicknessF = 3.0

   # Panel resources
   pnlres = Ngl.Resources()
   pnlres.nglFrame = False
   pnlres.nglPanelXWhiteSpacePercent = 5.0
   pnlres.nglPanelYWhiteSpacePercent = 5.0

   # Open workstation
   wks_type = "png"
   wks = Ngl.open_wks(wks_type, testfig + "/test_" + mpi + "-" + omp + "_" + suffix)

   # Plot
   plot = Ngl.xy(wks, x, mse_mean, xyres)

   # Add lines
   line = []
   for ifac in range(0, nfac):
      plres.gsMarkerColor = "black"
      line.append(Ngl.add_polyline(wks, plot, [x[ifac], x[ifac]], [mse_mean[0,ifac]-mse_std[0,ifac], mse_mean[0,ifac]+mse_std[0,ifac]], plres))
      plres.gsMarkerColor = "red"
      line.append(Ngl.add_polyline(wks, plot, [x[ifac], x[ifac]], [mse_mean[1,ifac]-mse_std[1,ifac], mse_mean[1,ifac]+mse_std[1,ifac]], plres))

   # Panel
   Ngl.draw(plot)

   # Advance frame
   Ngl.frame(wks)

   # Delete frame
   Ngl.delete_wks(wks)
