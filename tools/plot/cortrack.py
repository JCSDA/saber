#!/usr/bin/env python3

import argparse
import Ngl, Nio
import numpy as np
import os

def cortrack(testdata, test, mpi, omp, suffix):
   # Open file
   f = Nio.open_file(testdata + "/" + test + "/test_" + mpi + "-" + omp + "_" + suffix + ".nc")

   # Get lon/lat
   lon = f.variables["lon"][:]
   lat = f.variables["lat"][:]

   # Variables
   var_list = ["cor_00","cor_06","tracker_00_0","tracker_00_1","tracker_06_0"]
   suffix_list = ["","_tracker","_wind"]
   colors = {}
   colors[""] = "blue"
   colors["_tracker"] = "red"
   colors["_wind"] = "green"

   # Get number of levels
   nl0 = f.variables[var_list[0]][:,:].shape[0]

   # Get number of timeslots
   nts = f.variables["londir"][:].shape[0]

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
   cres.lbLabelBarOn = False
   cres.cnInfoLabelOn = False
   cres.cnLineLabelsOn = False
   cres.cnLinesOn = False
   cres.cnNoDataLabelOn = False
   cres.cnMissingValPerimOn = True
   cres.cnMissingValPerimColor = "black"
   cres.cnLevelSelectionMode = "ExplicitLevels"
   cres.cnLevels = Ngl.fspan(-1.0,1.0,21)
   cres.mpOutlineOn = False
   cres.mpLandFillColor = -1
   cres.mpProjection = "CylindricalEquidistant"
   cres.mpLimitMode = "LatLon"
   cres.mpMinLatF = 0.0
   cres.mpMaxLatF = 90.0
   cres.mpGridSpacingF = 45.0
   cres.sfXArray = lon
   cres.sfYArray = lat

   # Polyline resources
   plres = Ngl.Resources()
   plres.gsLineThicknessF = 3.0

   # Polymarker resources
   pmres = Ngl.Resources()
   pmres.gsMarkerIndex = 1
   pmres.gsMarkerSizeF = 0.03

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
      Ngl.define_colormap(wks, "BlWhRe")

      # Plots
      plot = []
      for il0 in range(0, nl0):
         plot.append(Ngl.contour_map(wks, field[il0,:], cres))

      # Panel
      Ngl.panel(wks, plot, [nl0,1], pnlres)

      # Advance frame
      Ngl.frame(wks)

      # Delete frame
      Ngl.delete_wks(wks)

   # Delete unnecessary resources
   del(cres.sfXArray)
   del(cres.sfYArray)
   del(cres.cnFillMode)
   del(cres.cnFillOn)
   del(cres.cnFillPalette)
   del(cres.cnInfoLabelOn)
   del(cres.cnLevelSelectionMode)
   del(cres.cnLevelSpacingF)
   del(cres.cnLineLabelsOn)
   del(cres.cnLinesOn)
   del(cres.cnMaxLevelValF)
   del(cres.cnMinLevelValF)
   del(cres.cnMissingValPerimColor)
   del(cres.cnMissingValPerimOn)
   del(cres.cnMonoFillColor)
   del(cres.cnMonoFillPattern)
   del(cres.cnNoDataLabelOn)
   del(cres.lbLabelBarOn)

   # Open workstation
   wks_type = "png"
   wks = Ngl.open_wks(wks_type, testfig + "/test_" + mpi + "-" + omp + "_" + suffix + "_coords")

   # Plot map
   map = Ngl.map(wks, cres)

   # Add lines and markers
   marker = []
   line = []
   for suffix in suffix_list:
      # Read lon/lat
      londir = f.variables["londir" + suffix][:]
      latdir = f.variables["latdir" + suffix][:]

      # Add lines
      plres.gsLineColor = colors[suffix]
      line.append(Ngl.add_polyline(wks, map, londir, latdir, plres))

      # Add markers
      for its in range(0, nts):
         if its==0:
            pmres.gsMarkerColor = "black"
         else:
            pmres.gsMarkerColor = colors[suffix]
         marker.append(Ngl.add_polymarker(wks, map, londir[its], latdir[its], pmres))

   # Draw map
   Ngl.draw(map)

   # Advance frame
   Ngl.frame(wks)

   # Delete frame
   Ngl.delete_wks(wks)
