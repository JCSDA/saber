#!/usr/bin/env python3

import argparse
import Ngl, Nio
import numpy as np
import numpy.ma as ma
import os

def sampling_grids(testdata, test, mpi, omp, suffix):
   # c3 plot
   ic3_plot = [0,2,4,8]

   # c2 plot
   ic2_plot = [0,2,4,8]

   # Loop over files
   first = True
   for impi in range(0, int(mpi)):
      # Open file
      f = Nio.open_file(testdata + "/" + test + "/test_" + mpi + "-" + omp + "_" + suffix + "_" + mpi.zfill(6) + "-" + format(impi+1, '06d') + ".nc")

      # Check what is in the file
      if first:
         sc3 = ("lon" in f.variables)
         local = ("lon_local" in f.variables)
         vbal = ("lon_vbal" in f.variables)

      # Get data
      if sc3:
         _FillValue = f.variables["lon"].attributes["_FillValue"]
      if local:
         _FillValue = f.variables["lon_local"].attributes["_FillValue"]
      if vbal:
         _FillValue = f.variables["lon_vbal"].attributes["_FillValue"]
      nc0a = f.variables["lon_c0a"].shape[0]
      if sc3:
         nc1a = f.variables["lon"].shape[2]
      if local:
         nc2a = f.variables["lon_local"].shape[1]
      if vbal:
         nc2b = f.variables["lon_vbal"].shape[1]
      if first:
         nl0 = f.variables["gmask_c0a"].shape[0]
         lon_c0 = np.zeros((0))
         lat_c0 = np.zeros((0))
         gmask_c0 = np.zeros((0,nl0))
         if sc3:
            nc3 = f.variables["lon"].shape[1]
            lon = np.zeros((0, nl0, nc3))
            lat = np.zeros((0, nl0, nc3))
         if local:
            nc1max_local = f.variables["lon_local"].shape[2]
            lon_local = np.zeros((0, nl0, nc1max_local))
            lat_local = np.zeros((0, nl0, nc1max_local))
         if vbal:
            nc1max_vbal = f.variables["lon_vbal"].shape[2]
            lon_vbal = np.zeros((0, nl0, nc1max_vbal))
            lat_vbal = np.zeros((0, nl0, nc1max_vbal))
         first = False
      for ic0a in range(0, nc0a):
         lon_c0 = np.append(lon_c0, f.variables["lon_c0a"][ic0a].data)
         lat_c0 = np.append(lat_c0, f.variables["lat_c0a"][ic0a].data)
         gmask_c0 = np.append(gmask_c0, [f.variables["gmask_c0a"][:,ic0a].data], axis=0)
      if sc3:
         for ic1a in range(0, nc1a):
            lon = np.append(lon, [f.variables["lon"][:,:,ic1a].data], axis=0)
            lat = np.append(lat, [f.variables["lat"][:,:,ic1a].data], axis=0)
      if local:
         for ic2a in range(0, nc2a):
            lon_local = np.append(lon_local, [f.variables["lon_local"][:,ic2a,:].data], axis=0)
            lat_local = np.append(lat_local, [f.variables["lat_local"][:,ic2a,:].data], axis=0)
      if vbal:
         for ic2b in range(0, nc2b):
            lon_vbal = np.append(lon_vbal, [f.variables["lon_vbal"][:,ic2b,:].data], axis=0)
            lat_vbal = np.append(lat_vbal, [f.variables["lat_vbal"][:,ic2b,:].data], axis=0)

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
   levs = f.attributes["nam_levs"].split(":")

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
   cres.sfXArray = lon_c0
   cres.sfYArray = lat_c0

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

   # Make output directory
   testfig = testdata + "/" + test + "/fig"
   if not os.path.exists(testfig):
      os.mkdir(testfig)

   if sc3:
      for ic3 in ic3_plot:
         # Open workstation
         wks_type = "png"
         wks = Ngl.open_wks(wks_type, testfig + "/test_" + mpi + "-" + omp + "_" + suffix + "_c3_" + str(ic3))

         # Plots
         plot = []
         for il0 in range(0, nl0):
            # Blank plot
            p = Ngl.contour_map(wks, gmask_c0[:,il0], cres)

            # Add polylines/polymarkers
            for ic1 in range(0, nc1):
                if (not lon.mask[ic1,il0,ic3]) and (not lat.mask[ic1,il0,ic3]):
                   Ngl.add_polyline(wks, p, [lon[ic1,il0,0], lon[ic1,il0,ic3]], [lat[ic1,il0,0], lat[ic1,il0,ic3]], plres)
                   Ngl.add_polymarker(wks, p, [lon[ic1,il0,0], lon[ic1,il0,ic3]], [lat[ic1,il0,0], lat[ic1,il0,ic3]], pmres)

            # Append plot
            plot.append(p)

         # Panel
         Ngl.panel(wks, plot, [nl0,1], pnlres)

         # Advance frame
         Ngl.frame(wks)

         # Delete frame
         Ngl.delete_wks(wks)

   if sc3 and (local or vbal):
      # Local sampling

      # Open workstation
      wks_type = "png"
      wks = Ngl.open_wks(wks_type, testfig + "/test_" + mpi + "-" + omp + "_" + suffix + "_c2_in_c1")

      # Plots
      plot = []
      for il0 in range(0, nl0):
         # Blank plot
         p = Ngl.contour_map(wks, gmask_c0[:,il0], cres)

         # Add polylines/polymarkers
         for ic1 in range(0, nc1):
            if (not lon.mask[ic1,il0,0]) and (not lat.mask[ic1,il0,0]):
                pmres.gsMarkerColor = "black"
                Ngl.add_polymarker(wks, p, lon[ic1,il0,0], lat[ic1,il0,0], pmres)
         if local:
            for ic2 in range(0, nc2):
               if (not lon_local.mask[ic2,il0,0]) and (not lat_local.mask[ic2,il0,0]):
                   pmres.gsMarkerColor = "red"
                   Ngl.add_polymarker(wks, p, lon_local[ic2,il0,0], lat_local[ic2,il0,0], pmres)
         if vbal:
            for ic2 in range(0, nc2btot):
               if (not lon_vbal.mask[ic2,il0,0]) and (not lat_vbal.mask[ic2,il0,0]):
                   pmres.gsMarkerColor = "red"
                   Ngl.add_polymarker(wks, p, lon_vbal[ic2,il0,0], lat_vbal[ic2,il0,0], pmres)

         # Append plot
         plot.append(p)

      # Panel
      Ngl.panel(wks, plot, [nl0,1], pnlres)

      # Advance frame
      Ngl.frame(wks)

      # Delete frame
      Ngl.delete_wks(wks)

      # Reset resource
      pmres.gsMarkerColor = "black"

   if local:
      # Local mask
      for ic2 in ic2_plot:
         # Open workstation
         wks_type = "png"
         wks = Ngl.open_wks(wks_type, testfig + "/test_" + mpi + "-" + omp + "_" + suffix + "_c2_" + str(ic2))

         # Plots
         plot = []
         for il0 in range(0, nl0):
            # Blank plot
            p = Ngl.contour_map(wks, gmask_c0[:,il0], cres)

            # Add polylines/polymarkers
            for ic1 in range(0, nc1max_local):
               if (not lon_local.mask[ic2,il0,ic1]) and (not lat_local.mask[ic2,il0,ic1]):
                  Ngl.add_polyline(wks, p, [lon_local[ic2,il0,0], lon_local[ic2,il0,ic1]], [lat_local[ic2,il0,0], lat_local[ic2,il0,ic1]], plres)
                  Ngl.add_polymarker(wks, p, [lon_local[ic2,il0,0], lon_local[ic2,il0,ic1]], [lat_local[ic2,il0,0], lat_local[ic2,il0,ic1]], pmres)

            # Append plot
            plot.append(p)

         # Panel
         Ngl.panel(wks, plot, [nl0,1], pnlres)

         # Advance frame
         Ngl.frame(wks)

         # Delete frame
         Ngl.delete_wks(wks)

   if vbal:
      # Reset resource
      pmres.gsMarkerColor = "black"

      # Vertical balance mask
      for ic2 in ic2_plot:
         # Open workstation
         wks_type = "png"
         wks = Ngl.open_wks(wks_type, testfig + "/test_" + mpi + "-" + omp + "_" + suffix + "_c2_" + str(ic2))

         # Plots
         plot = []
         for il0 in range(0, nl0):
            # Blank plot
            p = Ngl.contour_map(wks, gmask_c0[:,il0], cres)

            # Add polylines/polymarkers
            for ic1 in range(0, nc1max_vbal):
               if (not lon_vbal.mask[ic2,il0,ic1]) and (not lat_vbal.mask[ic2,il0,ic1]):
                  Ngl.add_polyline(wks, p, [lon_vbal[ic2,il0,0], lon_vbal[ic2,il0,ic1]], [lat_vbal[ic2,il0,0], lat_vbal[ic2,il0,ic1]], plres)
                  Ngl.add_polymarker(wks, p, [lon_vbal[ic2,il0,0], lon_vbal[ic2,il0,ic1]], [lat_vbal[ic2,il0,0], lat_vbal[ic2,il0,ic1]], pmres)

            # Append plot
            plot.append(p)

         # Panel
         Ngl.panel(wks, plot, [nl0,1], pnlres)

         # Advance frame
         Ngl.frame(wks)

         # Delete frame
         Ngl.delete_wks(wks)

