#!/usr/bin/env python3

import argparse
import Ngl, Nio
import numpy as np
import os

def umf(testdata, test, mpi, omp, suffix):
   # Open file
   f = Nio.open_file(testdata + "/" + test + "/test_" + mpi + "-" + omp + "_" + suffix + ".nc")

   # Get lon/lat
   lon = f.variables["lon"][:]
   lat = f.variables["lat"][:]

   # Variables
   var_list = ["m2","m4","kurt"]

   # Get number of levels
   nl0 = f.variables[var_list[0]][:,:].shape[0]

   # Contour resources
   cres = Ngl.Resources()
   cres.nglDraw = False
   cres.nglFrame = False
   cres.nglMaximize = True
   cres.cnFillOn = True
   cres.cnFillMode = "AreaFill"
   cres.trGridType = "TriangularMesh"
   cres.cnMonoFillPattern = True
   cres.cnMonoFillColor = False
   cres.cnFillPalette = "WhiteBlueGreenYellowRed"
   cres.lbLabelBarOn = False
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
   pnlres.nglPanelLabelBar = True

   # Make output directory
   testfig = testdata + "/" + test + "/fig"
   if not os.path.exists(testfig):
      os.mkdir(testfig)

   for var in var_list:
      # Read variable
      field = f.variables[var][:,:]

      # Open workstation
      wks_type = "png"
      wks = Ngl.open_wks(wks_type, testfig + "/test_" + mpi + "-" + omp + "_" + suffix + "_" + var)

      # Plots
      plot = []
      for il0 in range(0, nl0):
         plot.append(Ngl.contour_map(wks, field[il0,:], cres))

      # Panel
      Ngl.panel(wks, plot, [nl0,1], pnlres)

      # Advance frame
      Ngl.frame(wks)
