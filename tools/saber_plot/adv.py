#!/usr/bin/env python3

import argparse
import Ngl, Nio
import numpy as np
import numpy.ma as ma
import os

def adv(testdata, test, mpi, omp, suffix):
   # Open file
   f = Nio.open_file(testdata + "/" + test + "/test_" + mpi + "-" + omp + "_" + suffix + "_diag.nc")

   # Get variables
   lon_c2 = f.variables["lon_c2"][:,:]
   lat_c2 = f.variables["lat_c2"][:,:]
   lon_c2_raw = f.variables["lon_c2_raw"][:,:,:]
   lat_c2_raw = f.variables["lat_c2_raw"][:,:,:]
   dist_c2_raw = f.variables["dist_c2_raw"][:,:,:]
   valid_raw = f.variables["valid_raw"][:,:]
   dist_raw = f.variables["dist_raw"][:,:]
   lon_c2_flt = f.variables["lon_c2_flt"][:,:,:]
   lat_c2_flt = f.variables["lat_c2_flt"][:,:,:]
   dist_c2_flt = f.variables["dist_c2_flt"][:,:,:]
   valid_flt = f.variables["valid_flt"][:,:]
   dist_flt = f.variables["dist_flt"][:,:]
   rhflt = f.variables["rhflt"][:,:]
   score_loc = f.variables["score_loc"][:]
   score_adv = f.variables["score_adv"][:]

   # Get number of points
   nc2 = lon_c2_raw.shape[2]

   # Get number of levels
   nl0 = lon_c2_raw.shape[1]

   # Get number of timeslots
   nts = lon_c2_raw.shape[0]+1

   # Get levels
   levs = f.attributes["nam_levs"].split(":")

   # Get timeslots
   timeslots = f.attributes["nam_timeslots"].split(":")
   x = Ngl.fspan(1, nts, nts)

   # Open file
   f = Nio.open_file(testdata + "/" + test + "/test_" + mpi + "-" + omp + "_" + suffix + "_test.nc")

   # Get lon/lat
   lon = f.variables["lon"][:]
   lat = f.variables["lat"][:]

   # Variables
   root_list = []
   for var in f.variables:
      if var.find("_cor_loc")>0:
         root_list.append(var.replace("_cor_loc", ""))

   # Blank field
   blank = f.variables[root_list[0] + "_cor_loc"][:,:]

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
   cres.lbLabelBarOn = False
   cres.cnInfoLabelOn = False
   cres.cnLineLabelsOn = False
   cres.cnLinesOn = False
   cres.cnNoDataLabelOn = False
   cres.cnMissingValPerimOn = True
   cres.cnMissingValPerimColor = "black"
   cres.cnLevelSelectionMode = "ExplicitLevels"
   cres.cnLevels = [-1.0,2.0]
   cres.cnFillColors = ["white","white","white"]
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

   # XY resources
   xyres = Ngl.Resources()
   xyres.nglFrame = False
   xyres.nglDraw = False
   xyres.xyMarkLineMode = "Lines"
   xyres.xyLineColors = ["red", "blue"]
   xyres.xyLineThicknessF = 4.0
   xyres.trXMinF = 1
   xyres.trXMaxF = nts
   xyres.tmXBMode = "Explicit"
   xyres.tmXBValues = x
   xyres.tmXBLabels = timeslots
   xyres.tiXAxisString = "Timeslot"

   # Panel resources
   pnlres = Ngl.Resources()
   pnlres.nglFrame = False

   # Make output directory
   testfig = testdata + "/" + test + "/fig"
   if not os.path.exists(testfig):
      os.mkdir(testfig)

   # Max distance
   dist_max = np.max(dist_c2_raw)

   for its in range(0, nts):
      # Open workstation
      wks_type = "png"
      wks = Ngl.open_wks(wks_type, testfig + "/test_" + mpi + "-" + omp + "_" + suffix + "_raw_" + timeslots[its])
      Ngl.define_colormap(wks, "MPL_brg")

      # Plots
      plot = []
      for il0 in range(0, nl0):
         # Blank plot
         p = Ngl.contour_map(wks, blank[il0,:], cres)

         # Add polylines/polymarkers
         if its==0:
            for ic2 in range(0, nc2):
               pmres.gsMarkerColor = 2
               Ngl.add_polymarker(wks, p, lon_c2[il0,ic2], lat_c2[il0,ic2], pmres)
         else:
            for ic2 in range(0, nc2):
               if not dist_c2_raw.mask[its-1,il0,ic2]:
                  dist_color = np.min([2+int(dist_c2_raw[its-1,il0,ic2]/dist_max*(61.0-2.0)), 61])
                  plres.gsLineColor = dist_color
                  pmres.gsMarkerColor = dist_color
                  Ngl.add_polyline(wks, p, [lon_c2[il0,ic2], lon_c2_raw[its-1,il0,ic2]],[lat_c2[il0,ic2],lat_c2_raw[its-1,il0,ic2]], plres)
                  Ngl.add_polymarker(wks, p, lon_c2[il0,ic2], lat_c2[il0,ic2], pmres)

         # Append plot
         plot.append(p)

      # Panel
      Ngl.panel(wks, plot, [nl0,1], pnlres)

      # Advance frame
      Ngl.frame(wks)

      # Delete frame
      Ngl.delete_wks(wks)

   # Max distance
   dist_max = np.max(dist_c2_flt)

   for its in range(0, nts):
      # Open workstation
      wks_type = "png"
      wks = Ngl.open_wks(wks_type, testfig + "/test_" + mpi + "-" + omp + "_" + suffix + "_flt_" + timeslots[its])
      Ngl.define_colormap(wks, "MPL_brg")

      # Plots
      plot = []
      for il0 in range(0, nl0):
         # Blank plot
         p = Ngl.contour_map(wks, blank[il0,:], cres)

         # Add polylines/polymarkers
         if its==0:
            for ic2 in range(0, nc2):
               pmres.gsMarkerColor = 2
               Ngl.add_polymarker(wks, p, lon_c2[il0,ic2], lat_c2[il0,ic2], pmres)
         else:
            for ic2 in range(0, nc2):
               if not dist_c2_flt.mask[its-1,il0,ic2]:
                  dist_color = np.min([2+int(dist_c2_flt[its-1,il0,ic2]/dist_max*(61.0-2.0)), 61])
                  plres.gsLineColor = dist_color
                  pmres.gsMarkerColor = dist_color
                  Ngl.add_polyline(wks, p, [lon_c2[il0,ic2],lon_c2_flt[its-1,il0,ic2]], [lat_c2[il0,ic2],lat_c2_flt[its-1,il0,ic2]], plres)
                  Ngl.add_polymarker(wks, p, lon_c2[il0,ic2], lat_c2[il0,ic2], pmres)

         # Append plot
         plot.append(p)

      # Panel
      Ngl.panel(wks, plot, [nl0,1], pnlres)

      # Advance frame
      Ngl.frame(wks)

      # Delete frame
      Ngl.delete_wks(wks)

   # Reset resources
   del(cres.cnLevels)
   del(cres.cnFillColors)
   cres.cnLevels = cres.cnLevels = Ngl.fspan(-1.0,1.0,21)
   pnlres.nglPanelLabelBar = True

   for root in root_list:
      # Open workstation
      wks_type = "png"
      wks = Ngl.open_wks(wks_type, testfig + "/test_" + mpi + "-" + omp + "_" + suffix + "_" + root + "_cor_loc")
      Ngl.define_colormap(wks, "BlWhRe")

      # Plots
      plot = []
      for il0 in range(0, nl0):
         plot.append(Ngl.contour_map(wks, f.variables[root_list[0] + "_cor_loc"][il0,:], cres))

      # Panel
      Ngl.panel(wks, plot, [nl0,1], pnlres)

      # Advance frame
      Ngl.frame(wks)

      # Delete frame
      Ngl.delete_wks(wks)

      # Open workstation
      wks_type = "png"
      wks = Ngl.open_wks(wks_type, testfig + "/test_" + mpi + "-" + omp + "_" + suffix + "_" + root + "_cor_adv")
      Ngl.define_colormap(wks, "BlWhRe")

      # Plots
      plot = []
      for il0 in range(0, nl0):
         plot.append(Ngl.contour_map(wks, f.variables[root_list[0] + "_cor_adv"][il0,:], cres))

      # Panel
      Ngl.panel(wks, plot, [nl0,1], pnlres)

      # Advance frame
      Ngl.frame(wks)

      # Delete frame
      Ngl.delete_wks(wks)

   # Reset resources
   pnlres.nglPanelLabelBar = False
   pnlres.nglPanelXWhiteSpacePercent = 5.0
   pnlres.nglPanelYWhiteSpacePercent = 5.0

   # Open workstation
   wks_type = "png"
   wks = Ngl.open_wks(wks_type, testfig + "/test_" + mpi + "-" + omp + "_" + suffix)
   Ngl.define_colormap(wks, "MPL_brg")

   # Arrays
   valid = np.zeros((2, nts, nl0))
   dist = np.zeros((2, nts, nl0))
   rhflt_full = np.zeros((nts, nl0))
   valid[0,0,:] = 1.0
   valid[1,0,:] = 1.0
   valid[0,1:,:] = valid_flt
   valid[1,1:,:] = valid_raw
   dist[0,0,:] = 0.0
   dist[1,0,:] = 0.0
   dist[0,1:,:] = dist_flt
   dist[1,1:,:] = dist_raw
   rhflt_full[0,:] = 0.0
   rhflt_full[1:,:] = rhflt

   # Plots
   plot = []
   for il0 in range(0, nl0):
      # Valid points
      xyres.tiMainString = "Valid points at level " + levs[il0]
      xyres.trYMinF = 0.0
      xyres.trYMaxF = 1.1
      plot.append(Ngl.xy(wks, x, valid[:,:,il0], xyres))

      # Distance
      xyres.tiMainString = "Distance at level " + levs[il0]
      xyres.trYMinF = 0.0
      xyres.trYMaxF = np.max(dist)*1.1
      plot.append(Ngl.xy(wks, x, dist[:,:,il0], xyres))

      # Filtering
      xyres.tiMainString = "Filtering at level " + levs[il0]
      xyres.trYMinF = 0.0
      xyres.trYMaxF = np.max(rhflt_full)*1.1
      plot.append(Ngl.xy(wks, x, rhflt_full[:,il0], xyres))

   # Panel
   Ngl.panel(wks, plot, [nl0,3], pnlres)

   # Advance frame
   Ngl.frame(wks)

   # Delete frame
   Ngl.delete_wks(wks)

   # Open workstation
   wks_type = "png"
   wks = Ngl.open_wks(wks_type, testfig + "/test_" + mpi + "-" + omp + "_" + suffix + "_score")
   Ngl.define_colormap(wks, "MPL_brg")

   # Array
   score = np.zeros((2, nts))
   score[0,0] = 1.0
   score[1,0] = 1.0
   score[0,1:] = score_adv
   score[1,1:] = score_loc

   # Plot
   xyres.tiMainString = "Correlation at zero separation score"
   xyres.trYMinF = 0.0
   xyres.trYMaxF = 1.1
   plot = Ngl.xy(wks, x, score, xyres)

   # Draw
   Ngl.draw(plot)

   # Advance frame
   Ngl.frame(wks)

   # Delete frame
   Ngl.delete_wks(wks)

