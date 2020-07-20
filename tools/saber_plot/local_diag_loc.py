#!/usr/bin/env python3

import argparse
import Ngl, Nio
import numpy as np
import numpy.ma as ma
import os

def local_diag_loc(testdata, test, mpi, omp, suffix):
   # Open file
   f = Nio.open_file(testdata + "/" + test + "/test_" + mpi + "-" + omp + "_" + suffix + ".nc")

   # Get lon/lat
   lon = f.variables["lon"][:]
   lat = f.variables["lat"][:]

   # Variables
   var_list = []
   for var in f.variables:
      if var.find("coef_ens")>0 or var.find("coef_sta")>0 or var.find("fit_rh")>0 or var.find("fit_rv")>0:
         var_list.append(var)

   # Get number of levels
   nl0 = f.variables[var_list[0]][:,:].shape[0]

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
      Ngl.define_colormap(wks, "WhiteBlueGreenYellowRed")

      # Plots
      plot = []
      _FillValue = f.variables[var].attributes["_FillValue"]
      for il0 in range(0, nl0):
         if (np.any(field[il0,:] != _FillValue)):
            plot.append(Ngl.contour_map(wks, field[il0,:], cres))

      # Panel
      Ngl.panel(wks, plot, [nl0,1], pnlres)

      # Advance frame
      Ngl.frame(wks)

      # Delete frame
      Ngl.delete_wks(wks)
