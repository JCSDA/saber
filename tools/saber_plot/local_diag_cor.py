#!/usr/bin/env python3

import argparse
from netCDF4 import Dataset
import Ngl
import numpy as np
import numpy.ma as ma
import os

def local_diag_cor(testdata, test, mpi, omp, suffix, testfig):
   # Open file
   f = Dataset(testdata + "/" + test + "/test_" + mpi + "-" + omp + "_" + suffix + ".nc", "r", format="NETCDF4")

   # Get _FillValue
   _FillValue = f.__dict__["_FillValue"]

   # Get lon/lat
   lon = f["lon"][:]
   lat = f["lat"][:]

   # Get vertical unit
   vunit = f["vunit"][:,:]

   # Get number of levels
   nl0 = vunit.shape[0]

   # Contour resources
   cres = Ngl.Resources()
   cres.nglDraw = False
   cres.nglFrame = False
   cres.nglMaximize = True
   cres.cnFillOn = True
   cres.cnFillMode = "AreaFill"
   cres.cnConstFEnableFill = True
   cres.trGridType = "TriangularMesh"
   cres.cnMonoFillPattern = True
   cres.cnMonoFillColor = False
   cres.lbLabelBarOn = True
   cres.lbOrientation = "horizontal"
   cres.lbLabelFontHeightF  = 0.008
   cres.cnInfoLabelOn = False
   cres.cnLineLabelsOn = False
   cres.cnLinesOn = False
   cres.cnNoDataLabelOn = False
   cres.cnMissingValPerimOn = True
   cres.cnMissingValPerimColor = "black"
   cres.mpOutlineOn = False
   cres.mpLandFillColor = -1
   cres.mpProjection = "CylindricalEquidistant"
   cres.mpLimitMode = "LatLon"
   cres.mpMinLatF = 0.0
   cres.mpMaxLatF = 90.0
   cres.mpGridSpacingF = 45.0
   cres.sfXArray = lon
   cres.sfYArray = lat

   # Panel resources
   pnlres = Ngl.Resources()
   pnlres.nglFrame = False

   for group in f.groups:
      for var in f.groups[group].variables:
         # Read variable
         field = f.groups[group][var][:,:]

         # Open workstation
         wks_type = "png"
         wks = Ngl.open_wks(wks_type, testfig + "/test_" + mpi + "-" + omp + "_" + suffix + "_" + group + "_" + var)
         Ngl.define_colormap(wks, "WhiteBlueGreenYellowRed")

         # Plots
         plot = []
         for il0 in range(0, nl0):
            if np.any(field[il0,:] != _FillValue):
               plot.append(Ngl.contour_map(wks, field[il0,:], cres))

         # Panel
         Ngl.panel(wks, plot, [nl0,1], pnlres)

         # Advance frame
         Ngl.frame(wks)

         # Delete frame
         Ngl.delete_wks(wks)
