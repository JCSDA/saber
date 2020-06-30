#!/usr/bin/env python3

import argparse
import Ngl, Nio
import numpy as np
import os

def corstats(testdata, test, mpi, omp, suffix):
   # Open file
   f = Nio.open_file(testdata + "/" + test + "/test_" + mpi + "-" + omp + "_" + suffix + ".nc")

   # Get correlation maximum statistics
   cor_max_avg = f.variables["cor_max_avg"][:,:,:]
   cor_max_std = f.variables["cor_max_std"][:,:,:]

   # Get number of levels
   nl0 = cor_max_avg.shape[2]

   # Get number of variables
   nv = cor_max_avg.shape[1]

   # Get number of timeslots
   nts = cor_max_avg.shape[0]

   # Get levels
   levs = f.attributes["nam_levs"].split(":")

   # Get variables
   variables = f.attributes["nam_variables"].split(":")

   # Get timeslots
   timeslots = f.attributes["nam_timeslots"].split(":")
   x = Ngl.fspan(1, nts, nts)

   # XY resources
   xyres = Ngl.Resources()
   xyres.nglFrame = False
   xyres.nglDraw = False
   xyres.xyMarkLineMode = "Markers"
   xyres.xyMarker = 1
   xyres.xyMarkerColors = "black"
   xyres.xyMarkerSizeF = 0.05
   xyres.tiXAxisString = "Timeslot"
   xyres.tiYAxisString = "Correlation maximum"
   xyres.nglYRefLine = 1.0
   xyres.trXMinF = 0.5
   xyres.trXMaxF = float(nts)+0.5
   xyres.tmXBMode = "Explicit"
   xyres.tmXBValues = x
   xyres.tmXBLabels = timeslots
   xyres.trYMinF = -0.1
   xyres.trYMaxF = 1.1
   xyres.tmYLMode = "Explicit"
   xyres.tmYLValues = Ngl.fspan(0.0, 1.0, 6)
   xyres.tmYLLabels = xyres.tmYLValues
   xyres.tmYLMinorOn = True
   xyres.tmYLMinorValues = Ngl.fspan(-0.1, 1.1, 25)

   # Polyline resources
   plres = Ngl.Resources()
   plres.gsLineThicknessF = 3.0

   # Panel resources
   pnlres = Ngl.Resources()
   pnlres.nglFrame = False
   pnlres.nglPanelXWhiteSpacePercent = 5.0
   pnlres.nglPanelYWhiteSpacePercent = 5.0

   # Make output directory
   testfig = testdata + "/" + test + "/fig"
   if not os.path.exists(testfig):
      os.mkdir(testfig)

   # Open workstation
   wks_type = "png"
   wks = Ngl.open_wks(wks_type, testfig + "/test_" + mpi + "-" + omp + "_" + suffix)

   # Plots
   plot = []
   line = []
   for il0 in range(0, nl0):
      for iv in range(0, nv):
         # Title
         xyres.tiMainString = "Variable " + variables[iv] + " at level " + levs[il0]

         # Plot
         p = Ngl.xy(wks, x, cor_max_avg[:,iv,il0], xyres)

         # Add lines
         for its in range(0, nts):
            line.append(Ngl.add_polyline(wks, p, [x[its], x[its]], [cor_max_avg[its,iv,il0]-cor_max_std[its,iv,il0], cor_max_avg[its,iv,il0]+cor_max_std[its,iv,il0]], plres))

         # Append plot
         plot.append(p)

   # Panel
   Ngl.panel(wks, plot, [nv,nts], pnlres)

   # Advance frame
   Ngl.frame(wks)

   # Delete frame
   Ngl.delete_wks(wks)
